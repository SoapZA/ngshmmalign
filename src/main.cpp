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

int no_threads;

int main(int argc, const char* argv[])
{
	std::cout.imbue(std::locale("en_US.UTF-8"));

	// parameters
	std::string profile_filename;
	std::string output_filename;
	std::string rejects_filename;

	std::vector<std::string> input_files;
	background_rates params;
	bool write_unpaired = false;
	bool exhaustive = false;
	bool verbose = false;
	bool differentiate_match_state = false;

	clip_mode read_clip_mode;
	uint32_t min_mapped_length;
	uint64_t random_seed;

	/* set up program options */
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
		(",E", "Use full-exhaustive search, avoiding indexed lookup")
		(",X", "Replace general aligned state 'M' with '=' (match) and 'X' (mismatch) in CIGAR")
		("seed,s", boost::program_options::value<uint64_t>(&random_seed)->default_value(0), "Value of seed for deterministic run. A value of 0 will pick a random seed from some non-deterministic entropy source")
		("soft", "Soft-clip reads. Clipped bases will still be in the sequence in the alignment")
		("HARD", "Extreme Hard-clip reads. Do not write hard-clip in CIGAR, as if the hard-clipped bases never existed. Mutually exclusive with previous option")
		("verbose,v", "Show progress indicator while aligning")
		("ambig,a", boost::program_options::value<double>(&params.low_frequency_cutoff)->default_value(0.05, "0.05"), "Minimum frequency for calling ambiguous base")
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

	/* 0.0) parse program options */
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
	exhaustive = global_options.count("-E");
	verbose = global_options.count("verbose");

	if (global_options.count("soft") && global_options.count("HARD"))
	{
		std::cerr << "ERROR: You cannot have both '--soft' and '--HARD' enabled.\n";
		return EXIT_FAILURE;
	}
	read_clip_mode = (global_options.count("soft") ? clip_mode::soft : (global_options.count("HARD") ? clip_mode::HARD : clip_mode::hard));
	differentiate_match_state = global_options.count("-X");

	input_files = global_options["input-files"].as<std::vector<std::string>>();

	// set OpenMP properties
	omp_set_num_threads(no_threads);
	if (random_seed)
	{
		// deterministic
		omp_set_schedule(omp_sched_static, 200);
	}
	else
	{
		// non-deterministic
		omp_set_schedule(omp_sched_dynamic, 200);
	}

	if (!boost::filesystem::exists(profile_filename))
	{
		std::cerr << "ERROR: Reference file '" << profile_filename << "' does not exist!\n";
		return EXIT_FAILURE;
	}

	/* 0.1) create HMM aligner object */
	auto ngs_aligner = single_end_aligner<int32_t>::create_aligner_instance(write_unpaired, input_files, min_mapped_length, argc, argv);

	/* 1) load reads */
	ngs_aligner->load_reads(input_files);

	/* 2) load parameters */
	ngs_aligner->load_parameters(profile_filename, params);

	/* 3) sort reads */
	ngs_aligner->sort_reads();

	/* 4) perform alignment */
	ngs_aligner->perform_alignment(read_clip_mode, random_seed, exhaustive, verbose, differentiate_match_state);

	/* 5) write alignment to output */
	ngs_aligner->write_alignment_to_file(output_filename, rejects_filename);

	return 0;
}