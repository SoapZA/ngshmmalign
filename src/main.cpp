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

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "aligner.hpp"
#include "dna_array.hpp"
#include "fasta.hpp"
#include "hmmalign.hpp"
#include "reference.hpp"

int num_threads;

clip_mode read_entry::m_clip = clip_mode::soft;
boost::string_ref read_entry::RNAME;

int main(int argc, const char* argv[])
{
	try
	{
		std::cout.imbue(std::locale("en_US.UTF-8"));
	}
	catch (std::runtime_error)
	{
		std::cerr << "Warning: en_US.UTF-8 could not be imbued, this is likely due to a missing locale on your system\n";
		std::cout.imbue(std::locale::classic());
	}

	// parameters
	std::string profile_filename;
	std::string output_filename;
	std::string rejects_filename;
	std::string consensus_name;

	background_rates params;

	int32_t min_required_mapped_bases;
	uint64_t random_seed;

	/* set up program options */
	// program options
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "Print this help");

	// configuration options
	boost::program_options::options_description config("Configuration");
	config.add_options()
		(",r", boost::program_options::value<decltype(profile_filename)>(&profile_filename), "File containing the profile/MSA of the reference")
		(",R", boost::program_options::value<decltype(profile_filename)>(&profile_filename), "File containing the profile/MSA of the reference. Will perform a comprehensive parameter estimation using MAFFT. Mutually exclusive with -r option")
		(",o", boost::program_options::value<decltype(output_filename)>(&output_filename)->default_value("aln.sam"), "Filename where alignment will be written to")
		("wrong,w", boost::program_options::value<decltype(rejects_filename)>(&rejects_filename)->default_value("/dev/null"), "Filename where alignment will be written that are filtered (too short, unpaired)")
		(",t", boost::program_options::value<decltype(num_threads)>(&num_threads)->default_value(std::thread::hardware_concurrency()), "Number of threads to use for alignment. Defaults to number of logical cores found")
		(",l", "Do not clean up MAFFT temporary MSA files")
		(",E", "Use full-exhaustive search, avoiding indexed lookup")
		(",X", "Replace general aligned state 'M' with '=' (match) and 'X' (mismatch) in CIGAR")
		(",N", boost::program_options::value<decltype(consensus_name)>(&consensus_name)->default_value("CONSENSUS"), "Name of consensus reference contig that will be created")
		(",U", "Loci with ambiguous bases get their emission probabilities according to their allele frequencies. In practice this is undesirable, as it leads to systematic accumulation of gaps in homopolymeric regions with SNVs")
		("seed,s", boost::program_options::value<decltype(random_seed)>(&random_seed)->default_value(42), "Value of seed for deterministic run. A value of 0 will pick a random seed from some non-deterministic entropy source")
		("hard", "Hard-clip reads. Clipped bases will NOT be in the sequence in the alignment")
		("HARD", "Extreme Hard-clip reads. Do not write hard-clip in CIGAR, as if the hard-clipped bases never existed. Mutually exclusive with previous option")
		(",v", "Show progress indicator while aligning")

		(",M", boost::program_options::value<decltype(min_required_mapped_bases)>(&min_required_mapped_bases)->default_value(std::numeric_limits<decltype(min_required_mapped_bases)>::max(), "L * 0.8"), "Minimum mapped length of read")
		(",a", boost::program_options::value<decltype(params.low_frequency_cutoff)>(&params.low_frequency_cutoff)->default_value(0.05, "0.05"), "Minimum frequency for calling ambiguous base")

		("error", boost::program_options::value<decltype(params.error_rate)>(&params.error_rate)->default_value(0.005, "0.005"), "Global substitution probability")
		("go", boost::program_options::value<decltype(params.gap_open)>(&params.gap_open)->default_value(1e-4, "1e-4"), "Gap open probability")
		("ge", boost::program_options::value<decltype(params.gap_extend)>(&params.gap_extend)->default_value(0.30, "0.30"), "Gap extend probability")

		("io", boost::program_options::value<decltype(params.insert_open)>(&params.insert_open)->default_value(5e-5, "5e-5"), "Insert open probability")
		("ie", boost::program_options::value<decltype(params.insert_extend)>(&params.insert_extend)->default_value(0.50, "0.50"), "Insert extend probability")

		("ep", boost::program_options::value<decltype(params.end_prob)>(&params.end_prob)->default_value(MAGIC_NUMBER, "1/L"), "Jump to end probability; usually 1/L, where L is the average length of the reads")

		("lco", boost::program_options::value<decltype(params.left_clip_open)>(&params.left_clip_open)->default_value(0.10, "0.10"), "Left clip open probability")
		("lce", boost::program_options::value<decltype(params.left_clip_extend)>(&params.left_clip_extend)->default_value(0.90, "0.90"), "Left clip extend probability")

		("rco", boost::program_options::value<decltype(params.right_clip_open)>(&params.right_clip_open)->default_value(MAGIC_NUMBER, "lco/L"), "Right clip open probability")
		("rce", boost::program_options::value<decltype(params.right_clip_extend)>(&params.right_clip_extend)->default_value(0.90, "0.90"), "Right clip extend probability");

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
			std::cout << visible << std::endl;
			exit(EXIT_SUCCESS);
		}

		boost::program_options::notify(global_options);
	}
	catch (boost::program_options::required_option& e)
	{
		if (e.get_option_name() == "--input-files")
		{
			std::cerr << "ERROR: You have provided no input files. ngshmmalign takes either 1 (single-end) or 2 (paired-end) input file(s)." << std::endl;
		}
		else
		{
			std::cerr << "ERROR: " << e.what() << std::endl;
		}
		exit(EXIT_FAILURE);
	}
	catch (boost::program_options::error& e)
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
		exit(EXIT_FAILURE);
	}

	/* 0.1) assign parameter values */
	const bool exhaustive = global_options.count("-E");
	const bool verbose = global_options.count("-v");
	const bool differentiate_match_state = global_options.count("-X");
	const bool ambig_bases_unequal_weight = global_options.count("-U");
	const bool perform_hmm_learning = global_options.count("-R");
	const bool keep_mafft_files = global_options.count("-l");

	switch (global_options.count("-r") + perform_hmm_learning)
	{
		case 0:
			std::cerr << "ERROR: You need to specify either '-r' or '-R' for the reference." << std::endl;
			exit(EXIT_FAILURE);
			break;

		case 2:
			std::cerr << "ERROR: You cannot specify both '-r' and '-R'." << std::endl;
			exit(EXIT_FAILURE);
			break;

		default:
			break;
	}

	const char* env_var = std::getenv("MAFFT_BIN");
	const std::string mafft(env_var ? env_var : "mafft");
	if (perform_hmm_learning)
	{
		// try and find MAFFT
		// 1. try via environmental variable
		// 2. then try via PATH
		// std::cout << "LINE: " << (mafft + " --help > /dev/null 2>&1") << std::endl;

		int result = std::system((mafft + " --help > /dev/null 2>&1").c_str());

		switch (WEXITSTATUS(result))
		{
			case 0:
			case 1:
				std::cout << "Found MAFFT via " << (env_var ? "MAFFT_BIN" : "PATH") << std::endl;
				break;

			case 127:
				std::cerr << "ERROR: Couldn't find MAFFT." << std::endl;
				exit(EXIT_FAILURE);

			default:
				std::cerr << "ERROR: Unknown exit code returned." << std::endl;
				exit(EXIT_FAILURE);
		}
	}

	if (global_options.count("hard") && global_options.count("HARD"))
	{
		std::cerr << "ERROR: You cannot have both '--hard' and '--HARD' enabled." << std::endl;
		exit(EXIT_FAILURE);
	}
	const clip_mode read_clip_mode = (global_options.count("hard") ? clip_mode::hard : (global_options.count("HARD") ? clip_mode::HARD : clip_mode::soft));

	if (min_required_mapped_bases <= 0)
	{
		std::cerr << "ERROR: -M has to be strictly positive." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<std::string> input_files(global_options["input-files"].as<std::vector<std::string>>());

// set OpenMP number of threads
#ifdef _OPENMP
	omp_set_num_threads(num_threads);
#endif

	if (!boost::filesystem::exists(profile_filename))
	{
		std::cerr << "ERROR: Reference file '" << profile_filename << "' does not exist!" << std::endl;
		exit(EXIT_FAILURE);
	}
	boost::filesystem::path data_root_full_path(output_filename);
	std::string data_root;
	if (data_root_full_path.has_parent_path())
	{
		// output_filename contains a parent_path, i.e. a/b/c/aln.sam
		data_root = data_root_full_path.parent_path().string() + "/";
	}
	else
	{
		// output_filename contains NO parent_path, i.e. aln.sam
		data_root = "./";
		output_filename = data_root + output_filename;
	}

	boost::filesystem::path rejects_full_path(rejects_filename);
	rejects_filename = (rejects_full_path.has_parent_path() ? std::string() : data_root) + rejects_filename;

	/* 0.2) create HMM aligner object */
	auto ngs_aligner = single_end_aligner<int32_t>::create_aligner_instance(input_files, min_required_mapped_bases, argc, argv);

	/* 1) load reads */
	ngs_aligner->load_reads(input_files);

	/* 2) load parameters */
	ngs_aligner->load_parameters(profile_filename, params, ambig_bases_unequal_weight);

	/* 3) perform parameter estimation */
	if (perform_hmm_learning)
	{
		ngs_aligner->estimate_parameters(data_root, mafft, params, random_seed, verbose, keep_mafft_files, ambig_bases_unequal_weight);
	}

	/* 4) perform alignment */
	ngs_aligner->perform_alignment(exhaustive, verbose);

	/* 5) perform post-alignment processing */
	ngs_aligner->post_alignment_processing(differentiate_match_state, random_seed, params.low_frequency_cutoff, params.error_rate, ambig_bases_unequal_weight);

	/* 6) write alignment to output */
	ngs_aligner->write_alignment_to_file(read_clip_mode, consensus_name, data_root, output_filename, rejects_filename);

	return EXIT_SUCCESS;
}
