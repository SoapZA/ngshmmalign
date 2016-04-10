/*
 * Copyright (c) 2016 David Seifert
 * 	
 * This file is part of ngshmmalign
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <config.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <thread>
#include <tuple>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "dna_array.hpp"
#include "hmmalign.hpp"
#include "fasta.hpp"
#include "reference.hpp"
#include "aligner.hpp"

const int col_width = 8;
int no_threads;

int main(int argc, char* argv[])
{
	std::cout.imbue(std::locale("en_US.UTF-8"));

	// parameters
	std::string profile_filename;
	std::string output_filename;
	std::string rejects_filename;

	std::vector<std::string> input_files;
	background_rates<double> params;
	bool paired_end = false;
	bool write_unpaired = false;
	bool verbose = false;
	char clip_base;
	double ambig_threshold;
	uint32_t min_mapped_length;

	/* 1) set up program options */
	// program options
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "Print this help");

	// configuration options
	boost::program_options::options_description config("Configuration");
	config.add_options()
		("ref,r", boost::program_options::value<std::string>(&profile_filename)->required(), "File containing the profile/MSA of the reference")
		("out,o", boost::program_options::value<std::string>(&output_filename)->default_value("aln.sam"), "Filename where alignment will be written to")
		("wrong,w", boost::program_options::value<std::string>(&rejects_filename)->default_value("/dev/null"), "Filename where alignment will be written that are filtered (too short, unpaired)")
		("nthreads,t", boost::program_options::value<int32_t>(&no_threads)->default_value(std::thread::hardware_concurrency()), "Number of threads to use for alignment. Defaults to number of logical cores found")
		("unp", "Keep unpaired reads")
		("verbose,v", "Show progress indicator while aligning")
		("ambig,a", boost::program_options::value<double>(&ambig_threshold)->default_value(0.075, "0.075"), "Minimum frequency for calling ambiguous base")
		("minLen,M", boost::program_options::value<uint32_t>(&min_mapped_length)->default_value(MAGIC_NUMBER, "L * 0.8"), "Minimum mapped length of read")

		("error", boost::program_options::value<double>(&params.error_rate)->default_value(0.005, "0.005"), "Global substitution probability")
		("go", boost::program_options::value<double>(&params.gap_open)->default_value(1e-4, "1e-4"), "Gap open probability")
		("ge", boost::program_options::value<double>(&params.gap_extend)->default_value(0.30, "0.30"), "Gap extend probability")

		("io", boost::program_options::value<double>(&params.insert_open)->default_value(5e-5, "5e-5"), "Insert open probability")
		("ie", boost::program_options::value<double>(&params.insert_extend)->default_value(0.50, "0.50"), "Insert extend probability")

		("ep", boost::program_options::value<double>(&params.end_prob)->default_value(MAGIC_NUMBER, "1/L"), "Jump to end probability; usually 1/L, where L is the average length of the reads")

		("lco", boost::program_options::value<double>(&params.left_clip_open)->default_value(0.10, "0.10"), "Left clip open probability")
		("lce", boost::program_options::value<double>(&params.left_clip_extend)->default_value(0.90, "0.90"), "Left clip extend probability")

		("rco", boost::program_options::value<double>(&params.right_clip_open)->default_value(MAGIC_NUMBER, "lco/L"), "Right clip open probability")
		("rce", boost::program_options::value<double>(&params.right_clip_extend)->default_value(0.90, "0.90"), "Right clip extend probability");

	// hidden options, i.e., input files
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("input-files", boost::program_options::value<std::vector<std::string>>()->required(), "input files");

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(config).add(hidden);
	boost::program_options::options_description visible("Allowed options");
	visible.add(generic).add(config);
	boost::program_options::positional_options_description p;
	p.add("input-files", -1);
	boost::program_options::variables_map global_options;

	/* 2) parse program options */
	try
	{
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), global_options);

		// show help options
		if (global_options.count("help"))
		{
			std::cout << visible << "\n";
			return EXIT_SUCCESS;
		}

		boost::program_options::notify(global_options);
	}
	catch (boost::program_options::required_option& e)
	{
		if (e.get_option_name() == "--input-files")
			std::cerr << "ERROR: You have provided no input files. ngshmmalign takes either 1 (single-end) or 2 (paired-end) input file(s).\n";
		else
			std::cerr << "ERROR: " << e.what() << "\n";
		return EXIT_FAILURE;
	}
	catch (boost::program_options::error& e)
	{
		std::cerr << "ERROR: " << e.what() << "\n";
		return EXIT_FAILURE;
	}

	write_unpaired = global_options.count("unp");
	verbose = global_options.count("verbose");
	clip_base = (global_options.count("soft") ? 'S' : 'H');
	input_files = global_options["input-files"].as<std::vector<std::string>>();
	omp_set_num_threads(no_threads);

	switch (input_files.size())
	{
		case 1:
			paired_end = false;
			break;

		case 2:
			paired_end = true;
			break;

		default:
			std::cerr << "ERROR: You have provided too many input files. ngshmmalign takes either 1 (single-end) or 2 (paired-end) input file(s).\n";
			return EXIT_FAILURE;
			break;
	}

	/* 3) load references */
	if (!boost::filesystem::exists(profile_filename))
	{
		std::cerr << "ERROR: Reference file '" << profile_filename << "' does not exist!\n";
		return EXIT_FAILURE;
	}

	std::string ref_ext(boost::filesystem::extension(profile_filename));
	parameter_pack<double> hmm_parameters;
	std::tuple<std::vector<dna_array<double, 5>>, std::vector<double>, std::vector<double>> msa_profile;

	std::cout << "1) ";
	if ((ref_ext == ".fasta") || (ref_ext == ".fas") || (ref_ext == ".fna"))
	{
		std::cout << "Loading profile HMM parameters from FASTA file.\n";
		std::vector<reference_haplotype> references = fasta_read<reference_haplotype>(profile_filename);

		msa_profile = build_parameters_from_msa(references, ambig_threshold);
	}

	/* 4) create HMM object */
	std::unique_ptr<single_end_aligner<int32_t>> ngs_aligner((paired_end ? new paired_end_aligner<int32_t>(input_files[0], input_files[1], write_unpaired, min_mapped_length, argc, argv) : new single_end_aligner<int32_t>(input_files[0], min_mapped_length, argc, argv)));

	/* 5) sort reads */
	ngs_aligner->sort_reads();

	/* 6) set parameters */
	ngs_aligner->set_parameters(std::get<0>(msa_profile), std::get<1>(msa_profile), std::get<2>(msa_profile), params);

	/* 7) Perform alignments */
	ngs_aligner->perform_alignment(clip_base, verbose);

	/* 8) Write output */
	ngs_aligner->write_alignment_to_file(output_filename, rejects_filename);

	return 0;
}