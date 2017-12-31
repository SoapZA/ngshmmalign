#ifndef NGSHMMALIGN_ALIGNER_IMPL_HPP
#define NGSHMMALIGN_ALIGNER_IMPL_HPP

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

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <random>
#include <ratio>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/progress.hpp>

#include "aligner.hpp"
#include "debug.hpp"
#include "fastq.hpp"

namespace
{

fastq_entry::fastq_entry(
	const char* id_ptr,
	const std::size_t id_len,
	const char* seq_ptr,
	const std::size_t seq_len,
	const char* qual_ptr,
	const std::size_t qual_len,
	const bool second_in_pair)
	: m_id(id_ptr, id_len),
	  m_seq(seq_ptr, seq_len),
	  m_qual(qual_ptr, qual_len),
	  m_second_in_pair(second_in_pair)
{
}

read_entry::read_entry(
	const char* id_ptr,
	const std::size_t id_len,
	const char* seq_ptr,
	const std::size_t seq_len,
	const char* qual_ptr,
	const std::size_t qual_len,
	const bool second_in_pair)
	: m_fastq_record(
		  id_ptr,
		  id_len,
		  seq_ptr,
		  seq_len,
		  qual_ptr,
		  qual_len,
		  second_in_pair)
{
}

inline void read_entry::print_sequence_qual(std::ostream& output, bool clip_bases, boost::string_ref str) const noexcept
{
	const uint32_t left_clip = (clip_bases ? get_opt_aln().m_left_clip_length : 0);
	const uint32_t right_clip = (clip_bases ? get_opt_aln().m_right_clip_length : 0);

	std::ostream_iterator<decltype(str)::value_type> out_it(output);
	std::copy(str.cbegin() + left_clip, str.cend() - right_clip, out_it);
}

inline void read_entry::set_opt_aln(std::size_t i) noexcept
{
	assert(i < m_cooptimal_alignments.size());

	std::swap(m_cooptimal_alignments.front(), m_cooptimal_alignments[i]);
	m_cooptimal_alignments.erase(m_cooptimal_alignments.begin() + 1, m_cooptimal_alignments.end());

	assert(m_cooptimal_alignments.size() == 1);
}

inline minimal_alignment& read_entry::get_opt_aln() noexcept
{
	return const_cast<minimal_alignment&>(static_cast<const read_entry&>(*this).get_opt_aln());
}

inline const minimal_alignment& read_entry::get_opt_aln() const noexcept
{
	return m_cooptimal_alignments.front();
}

// output
std::ostream& operator<<(
	std::ostream& output,
	const read_entry& read) noexcept
{
	auto write_clip_CIGAR = [](std::ostream& output, std::size_t length, clip_mode clip) {
		switch (clip)
		{
			case clip_mode::soft:
				if (length)
				{
					output << length << 'S';
				}
				break;

			case clip_mode::hard:
				if (length)
				{
					output << length << 'H';
				}
				break;

			case clip_mode::HARD:
				break;
		}
	};

	///////////////////////////////////////
	// QNAME + FLAG + RNAME + POS + MAPQ //
	///////////////////////////////////////
	uint16_t FLAG = (read.m_paired_read * 0x1) | (read.m_mapped_in_proper_pair * 0x2) | (read.m_reverse_compl * 0x10) | (read.m_mate_reverse_compl * 0x20) | (read.m_paired_read * (read.m_fastq_record.m_second_in_pair ? 0x80 : 0x40));

	output
		<< read.m_fastq_record.m_id << '\t'
		<< FLAG << '\t'
		<< read.RNAME << '\t'
		<< read.get_opt_aln().m_POS + 1 << '\t'
		<< read.MAPQ << '\t';

	///////////
	// CIGAR //
	///////////
	write_clip_CIGAR(output, read.get_opt_aln().m_left_clip_length, read.m_clip);
	for (const auto& i : read.get_opt_aln().m_CIGAR)
	{
		output << i.second << i.first;
	}
	write_clip_CIGAR(output, read.get_opt_aln().m_right_clip_length, read.m_clip);

	//////////////////////////
	// RNEXT + PNEXT + TLEN //
	//////////////////////////
	output
		<< '\t' << read.RNEXT
		<< '\t' << read.PNEXT + 1
		<< '\t' << read.TLEN << '\t';

	////////////////
	// SEQ + QUAL //
	////////////////
	bool clip_bases = ((read.m_clip == clip_mode::hard) || (read.m_clip == clip_mode::HARD));
	read.print_sequence_qual(output, clip_bases, read.m_fastq_record.m_seq);
	output << '\t';
	read.print_sequence_qual(output, clip_bases, read.m_fastq_record.m_qual);

	/////////////////////////////////////
	// Alignment Score + Edit Distance //
	/////////////////////////////////////
	return output
		<< '\t' << "NM:i:" << read.NM
		<< '\t' << "MD:Z:" << read.MD_tag
		<< '\t' << "AS:i:" << read.SCORE;
}

template <typename T>
uint32_t single_end_aligner<T>::get_length_profile() const noexcept
{
	boost::accumulators::accumulator_set<int64_t, boost::accumulators::features<boost::accumulators::tag::median>> acc;
	for (std::vector<read_entry>::size_type i = 0; i < std::min(10000ul, m_reads.size()); ++i)
	{
		acc(m_reads[i].m_fastq_record.m_seq.length());
	}

	const int32_t median_length = boost::accumulators::median(acc);
	int32_t temp, diff = std::numeric_limits<int32_t>::max();
	uint32_t result = 25;

	for (auto i : { 25, 36, 50, 75, 150, 250, 300 }) // Illumina profiles
	{
		temp = std::abs(median_length - i);
		if (temp < diff)
		{
			diff = temp;
			result = i;
		}
	}

	return result;
}

// 1. factory function
template <typename T>
std::unique_ptr<single_end_aligner<T>> single_end_aligner<T>::create_aligner_instance(
	const std::vector<std::string>& input_files,
	const int32_t min_required_mapped_bases_,
	const int argc,
	const char** argv) noexcept
{
	switch (input_files.size())
	{
		case 0:
			std::cerr << "ERROR: You have provided no input files. ngshmmalign takes either 1 (single-end) or 2 (paired-end) input file(s)." << std::endl;
			exit(EXIT_FAILURE);
			break;

		case 1:
			return std::unique_ptr<single_end_aligner<T>>(new single_end_aligner<T>(min_required_mapped_bases_, argc, argv));
			break;

		case 2:
			return std::unique_ptr<single_end_aligner<T>>(new paired_end_aligner<T>(min_required_mapped_bases_, argc, argv));
			break;

		default:
			std::cerr << "ERROR: You have provided too many input files. ngshmmalign takes either 1 (single-end) or 2 (paired-end) input file(s)." << std::endl;
			exit(EXIT_FAILURE);
			break;
	}
}

template <typename T>
single_end_aligner<T>::single_end_aligner(
	const int32_t min_required_mapped_bases_,
	const int argc_,
	const char** argv_) noexcept
	: m_argc(argc_),
	  m_argv(argv_),
	  m_min_required_mapped_bases(min_required_mapped_bases_)
{
}

template <typename T>
paired_end_aligner<T>::paired_end_aligner(
	const int32_t min_required_mapped_bases_,
	const int argc_,
	const char** argv_) noexcept
	: single_end_aligner<T>(min_required_mapped_bases_, argc_, argv_)
{
}

// 2. load reads
template <typename T>
void single_end_aligner<T>::load_reads(
	const std::vector<std::string>& input_files) noexcept
{
	load_reads_impl(input_files);
}

template <typename T>
void single_end_aligner<T>::load_reads_impl(
	const std::vector<std::string>& input_files) noexcept
{
	// only a single file allowed for single-end reads
	assert(input_files.size() == 1);
	m_reads.clear();

	m_read_file_name = input_files.front();
	std::cout << ++m_phase << ") Loading single-end reads into memory from " << m_read_file_name << std::endl;
	fastq_read(m_read_file_name, m_reads, false, false);
}

template <typename T>
void paired_end_aligner<T>::load_reads_impl(
	const std::vector<std::string>& input_files) noexcept
{
	// need two files for paired-end reads
	assert(input_files.size() == 2);
	m_reads.clear();

	m_read_file_name = input_files.front();
	std::cout << ++m_phase << ") Loading (first) paired-end reads into memory from " << m_read_file_name << std::endl;
	fastq_read(m_read_file_name, m_reads, true, false);

	m_read_file_name2 = input_files.back();
	std::cout << "   Loading (second) paired-end reads into memory from " << m_read_file_name2 << std::endl;
	fastq_read(m_read_file_name2, m_reads, true, true);
}

// 3. load parameters
template <typename T>
void single_end_aligner<T>::load_parameters(
	const std::string& input_file,
	background_rates& error_rates,
	const bool ambig_bases_unequal_weight) noexcept
{
	std::cout << ++m_phase << ") Loading profile HMM parameters" << std::endl;

	uint32_t L = get_length_profile();
	std::cout << "   Using Illumina length profile " << L << std::endl
			  << "   Parameters for pHMM:" << std::endl;

	if (error_rates.end_prob == MAGIC_NUMBER)
	{
		error_rates.end_prob = 1.0 / L;
	}
	if (error_rates.right_clip_open == MAGIC_NUMBER)
	{
		error_rates.right_clip_open = error_rates.left_clip_open / L;
	}
	if (m_min_required_mapped_bases == std::numeric_limits<decltype(m_min_required_mapped_bases)>::max())
	{
		m_min_required_mapped_bases = (L * 4) / 5;
	}
	std::cout << std::endl
			  << error_rates << '\t' << "Minimum mapped length:   " << m_min_required_mapped_bases << std::endl
			  << std::endl;

	std::string ref_ext(boost::filesystem::extension(input_file));
	if ((ref_ext == ".fasta") || (ref_ext == ".fas") || (ref_ext == ".fna"))
	{
		std::cout << "   Input = FASTA" << std::endl;
		m_parameters.set_parameters(
			input_file,
			error_rates,
			L,
			ambig_bases_unequal_weight,
			typename reference_genome<T>::msa_input{});
	}
	else if (ref_ext == ".hmm")
	{
		std::cout << "   Input = " PACKAGE_NAME " serialised format" << std::endl;
		m_parameters.set_parameters(
			input_file,
			error_rates,
			L,
			ambig_bases_unequal_weight,
			typename reference_genome<T>::serialized_input{});
	}
	else
	{
		std::cerr << "ERROR: file ending" << ref_ext << " not recognized" << std::endl;
		exit(EXIT_FAILURE);
	}
}

// 4. perform parameter estimation
struct partioned_genome
{
	struct region
	{
		int32_t m_start;
		int32_t m_end;

		int32_t m_analysis_start;
		int32_t m_analysis_end;

		std::string m_files_dir;
		std::vector<const read_entry*> m_reads;

		region(const int32_t start_, const int32_t end_, const int32_t analysis_start_, const int32_t analysis_end_, const std::string& tmpdir, const int32_t num_digits) noexcept
			: m_start(start_),
			  m_end(end_),
			  m_analysis_start(analysis_start_),
			  m_analysis_end(analysis_end_)
		{
			std::ostringstream file_dir;
			file_dir
				<< tmpdir
				<< "/reads_"
				<< std::setw(num_digits) << std::setfill('0') << m_start
				<< "_"
				<< std::setw(num_digits) << std::setfill('0') << m_end;
			m_files_dir = file_dir.str();
		}
	};
	std::vector<region> m_regions;
	const int32_t m_offset;
	const int32_t m_num_bins;

	const int32_t m_min_required_mapped_bases;

	partioned_genome() = delete;
	partioned_genome(const partioned_genome&) = delete;
	partioned_genome(partioned_genome&&) = delete;
	partioned_genome& operator=(const partioned_genome& other) = delete;
	partioned_genome& operator=(partioned_genome&& other) = delete;

	partioned_genome(const std::string& tmpdir, const int32_t L, const int32_t read_profile) noexcept
		: m_offset(read_profile / 6),
		  m_num_bins(L / m_offset),
		  m_min_required_mapped_bases(read_profile * 0.8)
	{
		m_regions.reserve(m_num_bins);
		const int32_t num_digits = std::to_string(L).length();

		for (int32_t i = 0; i < L - read_profile; i += m_offset)
		{
			const int32_t end = std::min<int32_t>(L, i + m_offset + read_profile);
			m_regions.emplace_back(i, end, i + 3 * m_offset, i + 4 * m_offset, tmpdir, num_digits);
		}
	}

	void add_reads(const std::vector<read_entry>& reads, std::default_random_engine& generator) noexcept
	{
#ifndef NDEBUG
		uint64_t num_unique_assignments = 0;
		uint64_t num_multiple_assignments = 0;
		uint64_t num_no_assignments = 0;
#endif
		const int32_t max_index = m_regions.size() - 1;

		for (const auto& read : reads)
		{
			if (read.m_cooptimal_alignments.front().m_mapped_bases < m_min_required_mapped_bases)
			{
				// at least 80% of the profile length needs to be in mapped
				// bases in the read to be considered for the MSA.
				continue;
			}

			const auto start = read.m_cooptimal_alignments.front().m_POS - read.m_cooptimal_alignments.front().m_left_clip_length;
			const auto end = read.m_cooptimal_alignments.front().m_POS + read.m_cooptimal_alignments.front().m_segment_length + read.m_cooptimal_alignments.front().m_right_clip_length;

			std::vector<int32_t> candidate_regions;
			for (int32_t i = std::min<int32_t>(start / m_offset, max_index); i >= 0; --i)
			{
				if ((m_regions.at(i).m_start <= start) && (end <= m_regions.at(i).m_end))
				{
					candidate_regions.push_back(i);
				}
				else
				{
					break;
				}
			}

			if (candidate_regions.size())
			{
				// got some candidates
				int32_t rand_assignment = std::uniform_int_distribution<int32_t>(0, candidate_regions.size() - 1)(generator);
				m_regions[candidate_regions[rand_assignment]].m_reads.push_back(&read);

#ifndef NDEBUG
				if (candidate_regions.size() == 1)
				{
					++num_unique_assignments;
				}
				else
				{
					++num_multiple_assignments;
				}
			}
			else
			{
				++num_no_assignments;
#endif
			}
		}

#ifndef NDEBUG
		std::cerr
			<< "Unique:   " << num_unique_assignments << std::endl
			<< "Multiple: " << num_multiple_assignments << std::endl
			<< "None:     " << num_no_assignments << std::endl;
#endif
	}

	void shuffle_reads(const int32_t max_reads_to_sample, std::default_random_engine& generator) noexcept
	{
		for (auto& i : m_regions)
		{
			// 1. shuffle reads, for (sub)-sampling
			const int32_t num_reads = i.m_reads.size();
			if (num_reads > max_reads_to_sample)
			{
				for (int32_t j = 0; j < max_reads_to_sample; ++j)
				{
					const int32_t rand_swap = std::uniform_int_distribution<int32_t>(j, num_reads - 1)(generator);
					std::swap(i.m_reads[j], i.m_reads[rand_swap]);
				}
			}
		}
	}

	void write_out_regions(const std::string& reference, const int32_t max_reads_to_sample) const noexcept
	{
		for (const auto& i : m_regions)
		{
			boost::filesystem::create_directories(i.m_files_dir);

			const int32_t max_num_reads = std::min<int32_t>(i.m_reads.size(), max_reads_to_sample);
			std::ofstream output(i.m_files_dir + "/reads.fasta");
			for (int32_t j = 0; j < max_num_reads; ++j)
			{
				output << '>' << i.m_reads[j]->m_fastq_record.m_id << '\n';
				i.m_reads[j]->print_sequence_qual(output, true, i.m_reads[j]->m_fastq_record.m_seq);
				output << '\n';
			}
			output << ">REF\n"
				   << reference.substr(i.m_start, i.m_end - i.m_start) << '\n';
			output.close();
		}
	}

	void run_mafft_on_files(const std::string& mafft, const int32_t num_threads) const noexcept
	{
		std::cout << "Running MAFFT on region" << std::endl;
		for (const auto& i : m_regions)
		{
			std::ostringstream mafft_cmd_line;
			mafft_cmd_line
				<< mafft
				<< " --thread " << num_threads
				<< " --localpair --maxiterate 1000 --quiet "
				<< i.m_files_dir << "/reads.fasta > "
				<< i.m_files_dir << "/reads_aln.fasta";

			std::cout << std::right << std::setw(8) << i.m_start << " -> " << std::setw(6) << i.m_end << std::endl;
			int result = std::system(mafft_cmd_line.str().c_str());

			if (result != 0)
			{
				std::cerr << "ERROR: MAFFT exited with error code " << result << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	void construct_new_reference(
		const double min_base_cov,
		const double min_base_cutoff,
		std::vector<dna_array<double, 5>>& E_p,
		std::vector<double>& M_D_p,
		std::vector<double>& D_D_p) const noexcept
	{
		E_p.clear();
		E_p.reserve(1.2 * m_regions.back().m_end);

		M_D_p.clear();
		M_D_p.reserve(1.2 * m_regions.back().m_end);

		D_D_p.clear();
		D_D_p.reserve(1.2 * m_regions.back().m_end);

		//uint64_t r = 0;

		std::cout << "Processing MAFFT aligned files" << std::endl
				  << std::fixed;
		for (const auto& i : m_regions)
		{
			std::cout << std::right << std::setw(8) << i.m_start << " -> " << std::setw(6) << i.m_end << std::endl;

			// 0. load aligned MSA
			std::vector<reference_haplotype> msa_reads(fasta_read<reference_haplotype>(i.m_files_dir + "/reads_aln.fasta"));
			if (msa_reads.empty())
			{
				// no reads in alignment, have to load reference
				msa_reads = fasta_read<reference_haplotype>(i.m_files_dir + "/reads.fasta");
			}
			const reference_haplotype ref_seq(msa_reads.back());

			if (msa_reads.size() > 1)
			{
				msa_reads.pop_back();
			}

			if (ref_seq.name != "REF")
			{
				std::cerr << "ERROR: Last sequence in " << i.m_files_dir << "/reads_aln.fasta is not the reference!" << std::endl;
				exit(EXIT_FAILURE);
			}

			// 1. collect statistics
			const int32_t total_length_of_MSA = ref_seq.sequence.length();
			std::vector<uint32_t> coverage(total_length_of_MSA, 0);
			std::vector<uint32_t> nongap_histogram(total_length_of_MSA, 0);
			std::vector<dna_array<double, 5>> base_dist(total_length_of_MSA, { 0.0, 0.0, 0.0, 0.0, 0.0 });

			for (int32_t j = 0; j < total_length_of_MSA; ++j)
			{
				for (const auto& k : msa_reads)
				{
					if ((k.start <= j) && (j < k.end))
					{
						++coverage[j];
						char base = k.sequence[j];
						if (base != '-')
						{
							// base
							++nongap_histogram[j];

							const auto it = wobble_to_ambig_bases.find(base);
							if (it == wobble_to_ambig_bases.end())
							{
								std::cerr << "ERROR: Could not map base '" << base << "' to unambiguous bases!" << std::endl;
								exit(EXIT_FAILURE);
							}

							const std::string& all_bases = it->second;
							for (char non_ambig_base : all_bases)
							{
								base_dist[j][non_ambig_base] += (1.0 / all_bases.length());
							}
						}
					}
				}
			}

			// 2. determine whether a locus is of high coverage
			std::vector<uint8_t> is_locus_cov_high;
			is_locus_cov_high.reserve(total_length_of_MSA);
			for (int32_t j = 0; j < total_length_of_MSA; ++j)
			{
				is_locus_cov_high.push_back(nongap_histogram[j] > 0.5 * coverage[j]);
			}

			// 3. find start of region with respect to reference
			// The offsets form a CLOSED interval, which is unlike
			// what ngshmmalign usually uses, namely half-open intervals
			int32_t ref_start = 0;
			if (i.m_start != m_regions.front().m_start)
			{
				// not the first region
				int32_t found_bases = 0;
				while (found_bases <= (i.m_analysis_start - i.m_start))
				{
					found_bases += (ref_seq.sequence[ref_start] != '-');
					++ref_start;
				}
				--ref_start;

				// start index has to be on a valid (high coverage locus)
				// otherwise, bases might become lost
				while (is_locus_cov_high[ref_start] == false)
				{
					++ref_start;
				}
			}

			int32_t ref_end = total_length_of_MSA - 1;
			if (i.m_end != m_regions.back().m_end)
			{
				// not the last region
				int32_t found_bases = 0;
				while (found_bases <= (i.m_end - i.m_analysis_end - 1))
				{
					found_bases += (ref_seq.sequence[ref_end] != '-');
					--ref_end;
				}
				++ref_end;

				// end index has to be on a valid (high coverage locus)
				// otherwise, bases might become lost
				while (is_locus_cov_high[ref_end] == false)
				{
					++ref_end;
				}
			}

			// 4. find all valid loci in the MSA region
			// TODO: add more adaptive and intelligent
			// homopolymer detector
			std::vector<int32_t> valid_loci;
			valid_loci.reserve(ref_end - ref_start + 1);
			for (int32_t j = ref_start; j <= ref_end; ++j)
			{
				if (coverage[j])
				{
					// A general problem of our estimation scheme is that it is rather sensitive
					// to artefacts in low coverage regimes. In order to decide whether a locus
					// in the alignment is 'valid', we therefore require at least the 2 * reciprocal
					// of 'min_base_cov' number of reads supporting the locus AND a minimum
					// fraction of 'min_base_cov', and if we have less, 50% must support the position.
					if (coverage[j] >= 2.0 / min_base_cov)
					{
						if (nongap_histogram[j] >= min_base_cov * coverage[j])
						{
							valid_loci.push_back(j);
						}
					}
					else
					{
						if (nongap_histogram[j] > 0.5 * coverage[j])
						{
							valid_loci.push_back(j);
						}
					}
				}
			}

			// 5. determine parameters
			auto remove_low_frequency_bases = [min_base_cutoff](dna_array<double, 5>& allel_freq) -> void {
				const double coverage = std::accumulate(&allel_freq[static_cast<std::size_t>(0)], &allel_freq[static_cast<std::size_t>(4)], 0.0);
				for (char k : { 'A', 'C', 'G', 'T' })
				{
					if (allel_freq[k] < min_base_cutoff * coverage)
					{
						allel_freq[k] = 0;
					}
				}
			};

			const int32_t num_valid_loci = valid_loci.size() - 1;
			for (int32_t j = 0; j < num_valid_loci; ++j)
			{
				const int32_t left_locus = valid_loci[j];
				const int32_t right_locus = valid_loci[j + 1];

				double num_M = 0;
				double num_M_D = 0;

				double num_D = 0;
				double num_D_D = 0;

				for (const auto& k : msa_reads)
				{
					if ((k.start <= left_locus) && (right_locus < k.end))
					{
						double& count = ((k.sequence[left_locus] == '-') ? (++num_D, num_D_D) : (++num_M, num_M_D));
						count += (k.sequence[right_locus] == '-');
					}
				}

				num_M += (!num_M);
				if ((num_M_D < 0.05 * num_M) || (num_M < 10))
				{
					num_M_D = 0;
				}

				num_D += (!num_D);
				if ((num_D_D < 0.05 * num_D) || (num_D < 10))
				{
					num_D_D = 0;
				}

				remove_low_frequency_bases(base_dist[left_locus]);
				E_p.push_back(base_dist[left_locus]);
				M_D_p.push_back(num_M_D / num_M);
				D_D_p.push_back(num_D_D / num_D);

				if ((i.m_end == m_regions.back().m_end) && (right_locus == valid_loci.back()))
				{
					// deal with the last global position
					remove_low_frequency_bases(base_dist[right_locus]);
					E_p.push_back(base_dist[right_locus]);
				}
			}
		}
	}
};

template <typename T>
void single_end_aligner<T>::estimate_parameters(
	const std::string& data_root,
	const std::string& mafft,
	const background_rates& error_rates,
	const uint64_t seed,
	const bool verbose,
	const bool keep_mafft_files,
	const bool ambig_bases_unequal_weight,
	const int32_t num_threads,
	const int32_t chunk_size) noexcept
{
	// 0. create random TMP directory name in order to avoid clashes
	std::default_random_engine rng(seed);
	std::random_device rd;
	std::ostringstream tmp_root_constr;
	tmp_root_constr
		<< data_root
		<< "TMP_"
		<< (seed ? std::uniform_int_distribution<uint16_t>()(rng) : std::uniform_int_distribution<uint16_t>()(rd));
	const std::string tmpdir_root(tmp_root_constr.str());

	// 1. perform rough alignment first
	std::cout << ++m_phase << ") Performing first rough alignment" << std::endl;
	perform_alignment_impl(false, verbose, num_threads, chunk_size);

	// 2. chop up genome and extract reads
	partioned_genome separated_reads(tmpdir_root, m_parameters.m_L, m_parameters.read_length_profile);
	separated_reads.add_reads(m_reads, rng);

	// 3. shuffle reads for subsampling
	constexpr int32_t max_reads = 500;
	separated_reads.shuffle_reads(max_reads, rng);

	// 4. write reads to FASTA files
	separated_reads.write_out_regions(m_parameters.m_ambig_ref, max_reads);

	// 5. run MAFFT on regions
	separated_reads.run_mafft_on_files(mafft, num_threads);

	// 6. parse MAFFT alignments
	std::vector<dna_array<double, 5>> E_p;
	std::vector<double> M_D_p;
	std::vector<double> D_D_p;
	separated_reads.construct_new_reference(
		0.05,
		error_rates.low_frequency_cutoff,
		E_p,
		M_D_p,
		D_D_p);

	m_parameters.set_parameters(std::move(E_p), std::move(M_D_p), std::move(D_D_p), error_rates, ambig_bases_unequal_weight);

	// 7. (optional) cleanup temporary MAFFT files
	if (keep_mafft_files == false)
	{
		boost::filesystem::remove_all(tmpdir_root);
	}
}

// 5. perform alignment
void display_progress(
	const std::atomic<int32_t>& global_idx,
	const int32_t total) noexcept
{
	boost::progress_display show_progress(total);
	int32_t now{};
	int32_t previous{};

	do
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(250));
		show_progress += (now - previous);

		previous = now;
		now = global_idx.load(std::memory_order_relaxed);
	} while (now < total);
}

template <typename T>
void align_thread_worker(
	const bool exhaustive,
	std::atomic<int32_t>& global_idx,
	std::vector<read_entry>& reads,
	std::mutex& mutex_lock,
	const int32_t chunk_size,
	const reference_genome<T> parameters,
	const int32_t max_read_length)
{
	// minimum number of matched k-mers of a read to be considered the correct strand
	constexpr const uint32_t min_hash_pos_better = 20;

	// maximum number of matched k-mers of a read to be considered the correct strand
	constexpr const uint32_t max_hash_pos_worse = 2;

	const int32_t total = reads.size();
	int32_t thread_idx{};

	// 1. store alignment results in temporary vector
	//    in order to avoid false sharing
	struct temp_result
	{
		read_entry* read;
		bool reverse_complement;
		int64_t score;
		std::vector<minimal_alignment> alignments;

		temp_result(
			read_entry* read_,
			bool reverse_complement_,
			int64_t score_,
			std::vector<minimal_alignment>&& alignments_) : read{ read_ },
															reverse_complement{ reverse_complement_ },
															score{ score_ },
															alignments{ std::move(alignments_) } {}
	};
	std::vector<temp_result> temp_result_alignments;

	while ((thread_idx = global_idx.fetch_add(chunk_size, std::memory_order_relaxed)) < total)
	{
		hmmalign<T> aligner;

		const int32_t idx_end = std::min(total, thread_idx + chunk_size);
		for (int32_t i = thread_idx; i < idx_end; ++i)
		{
			const read_entry& read_entry_ref = reads[i];
			const boost::string_ref seq = read_entry_ref.m_fastq_record.m_seq;

			// 1. first use heuristics to determine strand
			T forward_score, reverse_score;
			bool align_forward = false;
			forward_score = std::numeric_limits<decltype(forward_score)>::min();
			auto forward_stats = parameters.find_pos(seq);

			bool align_reverse = false;
			reverse_score = std::numeric_limits<decltype(reverse_score)>::min();

			std::string rev_seq;
			rev_seq.resize(seq.size());
			std::transform(seq.crbegin(), seq.crend(), rev_seq.begin(), rev_comp_char);

			auto reverse_stats = parameters.find_pos(rev_seq);

			bool exhaustive_alignment = exhaustive;

			if ((forward_stats.num_samples > min_hash_pos_better) && (reverse_stats.num_samples < max_hash_pos_worse))
			{
				// probably forward
				align_forward = true;
			}
			else
			{
				if ((reverse_stats.num_samples > min_hash_pos_better) && (forward_stats.num_samples < max_hash_pos_worse))
				{
					// probably reverse
					align_reverse = true;
				}
				else
				{
					// couldn't determine strandness, have to do full exhaustive
					align_forward = true;
					align_reverse = true;

					exhaustive_alignment = true;
				}
			}

			// 2. perform alignments per strand
			auto per_strand_alignment = [
#ifndef NDEBUG
											&alignment_stats,
#endif
											&aligner,
											parameters](const boost::string_ref& query, bool exhaustive, const typename reference_genome<T>::index_stat& heuristics, std::vector<minimal_alignment>& alignments) -> T {
				int32_t POS;
				int32_t end_POS;
				T score;

				exhaustive |= (heuristics.num_samples <= min_hash_pos_better);
				exhaustive |= (heuristics.sd > 1000);

				if (exhaustive)
				{
					POS = 0;
					end_POS = parameters.m_L;
#ifndef NDEBUG
					++alignment_stats.num_exhaustive_alignments;
#endif
				}
				else
				{
					// heuristic buffer to add left and right of
					// supposed alignment region
					int32_t buffer = 100 + (6 * heuristics.sd) / 5;

					POS = heuristics.POS - buffer;
					POS = std::max<int32_t>(0, POS);

					end_POS = heuristics.POS + query.length() + buffer;
					end_POS = std::min<int32_t>(end_POS, parameters.m_L);

#ifndef NDEBUG
					// // test for detection of failed heuristic
					// POS = heuristics.POS + 100;
					// POS = std::max<int32_t>(0, POS);
					//
					// end_POS = heuristics.POS + query.length() - 100;
					// end_POS = std::min<int32_t>(end_POS, parameters.m_L);
					++alignment_stats.num_indexed_alignments;
#endif
				}

				alignments.clear();
				score = aligner.viterbi(parameters, query, POS, end_POS, alignments);

				DEBUG_TRACE(std::cerr
					<< "POS:" << std::right
					<< std::setw(6) << POS
					<< std::setw(6) << "End:"
					<< std::setw(6) << end_POS
					<< std::setw(14) << "True POS:"
					<< std::setw(6) << alignments.front().m_POS
					<< std::setw(14) << "True End:"
					<< std::setw(6) << alignments.front().m_POS + alignments.front().m_segment_length << std::endl);

				if ((((alignments.front().m_POS - alignments.front().m_left_clip_length < POS) && (POS != 0)) || ((alignments.front().m_POS + alignments.front().m_segment_length + alignments.front().m_right_clip_length > end_POS) && (end_POS != static_cast<int32_t>(parameters.m_L))))
					&& (exhaustive == false))
				{
					// window size probably too small, perform sacrificial full-length alignment
					alignments.clear();
					score = aligner.viterbi(parameters, query, 0, parameters.m_L, alignments);

#ifndef NDEBUG
					++alignment_stats.num_sacrificial_alignments;
					std::cerr << "Had to perform sacrificial exhaustive alignment" << std::endl;
#endif
				}

				return score;
			};

			std::vector<minimal_alignment> forward_alns, reverse_alns;
			if (align_forward)
			{
				forward_score = per_strand_alignment(seq, exhaustive_alignment, forward_stats, forward_alns);
			}

			if (align_reverse)
			{
				reverse_score = per_strand_alignment(rev_seq, exhaustive_alignment, reverse_stats, reverse_alns);
			}

			// 3. pick correct alignment
			if (reverse_score > forward_score)
			{
				// reverse complementary is better
				temp_result_alignments.emplace_back(&reads[i], true, reverse_score, std::move(reverse_alns));
			}
			else
			{
				// forward is better
				temp_result_alignments.emplace_back(&reads[i], false, forward_score, std::move(forward_alns));
			}
		}
	}

	// 2. write temporary results back into main vector
	std::lock_guard<std::mutex> lock{ mutex_lock };
	for (auto& j : temp_result_alignments)
	{
		read_entry& read_entry_ref = *(j.read);

		read_entry_ref.SCORE = j.score;
		read_entry_ref.m_cooptimal_alignments = std::move(j.alignments);

		if (j.reverse_complement)
		{
			// flip reverse complementary bit
			// setting it to true is not enough
			// imagine:
			// 1st run:
			//     <-------    Sequence + QUAL get reverse complemented
			//
			// 2nd run:
			//     ------->    Sequence + QUAL are forward aligned
			//                 but have originally been reverse complemented
			read_entry_ref.m_reverse_compl = !(read_entry_ref.m_reverse_compl);

			std::string& seq = read_entry_ref.m_fastq_record.m_seq;
			std::string rev_seq;
			rev_seq.resize(seq.size());
			std::transform(seq.crbegin(), seq.crend(), rev_seq.begin(), rev_comp_char);
			seq = std::move(rev_seq);

			std::reverse(read_entry_ref.m_fastq_record.m_qual.begin(), read_entry_ref.m_fastq_record.m_qual.end());
		}
	}
};

template <typename T>
void single_end_aligner<T>::perform_alignment_impl(
	const bool exhaustive,
	const bool verbose,
	const int32_t num_threads,
	const int32_t chunk_size) noexcept
{
	static_assert(std::is_integral<T>::value, "T needs to be an integral type!\n");

	std::cout << "   Aligning reads (" << m_reads.size() << ")" << std::flush;

#ifndef NDEBUG
	struct alignment_stats_struct
	{
		uint64_t num_exhaustive_alignments = 0;
		uint64_t num_indexed_alignments = 0;
		uint64_t num_sacrificial_alignments = 0;
	} alignment_stats;
	std::cerr << std::endl;
#endif

	std::vector<std::thread> aln_threads;
	aln_threads.reserve(num_threads + 1);

	std::atomic<int32_t> global_idx{ 0 };
	std::mutex global_lock;

	if (verbose)
	{
		aln_threads.emplace_back(display_progress, std::cref(global_idx), m_reads.size());
	}

	for (int32_t i = 0; i < num_threads; ++i)
	{
		aln_threads.emplace_back(std::thread(align_thread_worker<T>, exhaustive, std::ref(global_idx), std::ref(m_reads), std::ref(global_lock), chunk_size, m_parameters, m_max_read_length.value()));
	}

	// wait for threads to complete
	for (auto& j : aln_threads)
	{
		j.join();
	}

#ifndef NDEBUG
	std::cerr
		<< "Number of exhaustive alignments:  " << alignment_stats.num_exhaustive_alignments << std::endl
		<< "Number of indexed alignments:     " << alignment_stats.num_indexed_alignments << std::endl
		<< "Number of sacrificial alignments: " << alignment_stats.num_sacrificial_alignments << std::endl;
#endif

	std::cout << std::endl;
}

template <typename T>
void single_end_aligner<T>::perform_alignment(
	const bool exhaustive,
	const bool verbose,
	const int32_t num_threads,
	const int32_t chunk_size) noexcept
{
	static_assert(std::is_integral<T>::value, "T needs to be an integral type!\n");

	std::cout << ++m_phase << ") Performing alignment" << std::endl
			  << "   Number of threads -t:    " << num_threads << std::endl;
	std::chrono::steady_clock::time_point time_start = std::chrono::steady_clock::now();

	// perform actual alignment
	perform_alignment_impl(exhaustive, verbose, num_threads, chunk_size);

	std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
	const double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_start).count() / 1E9;
	const double performance = m_reads.size() / duration;

	std::cout << std::endl
			  << "   Alignment took "
			  << static_cast<uint32_t>(duration) / 3600 << "h "
			  << (static_cast<uint32_t>(duration) / 60) % 60 << "min "
			  << static_cast<uint32_t>(duration) % 60 << "s ("
			  << std::fixed << std::setprecision(1)
			  << performance / num_threads << " reads/thread/s, " << performance << " reads/s)." << std::endl;
}

// 6. post-alignment processing
template <typename T>
void single_end_aligner<T>::post_alignment_processing(const bool differentiate_match_state, const uint64_t seed, const double min_freq, const double error_rate, bool ambig_bases_unequal_weight) noexcept
{
	std::cout << ++m_phase << ") Performing post-alignment processing" << std::endl;

	// a. check whether all co-optimal alignments have at least m_min_required_mapped_bases bases
	std::cout << "   Checking for reads with too few mapped bases" << std::endl;

	int32_t min_mapped_bases;
	for (auto& i : m_reads)
	{
		i.m_invalid = false;
		min_mapped_bases = std::numeric_limits<decltype(min_mapped_bases)>::max();

		// determine minimum mapped bases for all co-optimal alignments
		for (const auto& j : i.m_cooptimal_alignments)
		{
			min_mapped_bases = std::min<decltype(min_mapped_bases)>(min_mapped_bases, j.m_mapped_bases);

			if (min_mapped_bases < m_min_required_mapped_bases)
			{
				// at least one alignment has less than m_min_required_mapped_bases mapped bases
				i.m_invalid = true;
				break;
			}
		}
	}

	// b. distribute gaps (and write SAM entries)
	std::cout << "   Distributing gaps among co-optimal alignments" << std::endl;

	std::vector<std::size_t> gap_histogram(m_parameters.m_L, 0);
	for (auto& i : m_reads)
	{
		int32_t pos_on_genome;
		std::size_t cur_gap_aln_score, best_gap_aln = 0;
		std::size_t best_gap_aln_score = std::numeric_limits<decltype(best_gap_aln_score)>::max();

		// 1. determine best alignment with respect to gaps
		for (std::size_t j = 0; j < i.m_cooptimal_alignments.size(); ++j)
		{
			cur_gap_aln_score = 0;
			pos_on_genome = i.m_cooptimal_alignments[j].m_POS;

			for (const auto& k : i.m_cooptimal_alignments[j].m_CIGAR)
			{
				switch (k.first)
				{
					case 'D':
						for (std::size_t l = 0; l < k.second; ++l)
						{
							cur_gap_aln_score += (gap_histogram[pos_on_genome + l]);
						}
					case 'M':
						pos_on_genome += k.second;

					default:
						break;
				}
			}

			if (cur_gap_aln_score < best_gap_aln_score)
			{
				best_gap_aln = j;
				best_gap_aln_score = cur_gap_aln_score;
			}
		}

		// 2. choose best alignment
		i.set_opt_aln(best_gap_aln);

		// 3. increment gap histogram
		if (i.m_invalid == false)
		{
			pos_on_genome = i.get_opt_aln().m_POS;
			for (const auto& k : i.get_opt_aln().m_CIGAR)
			{
				switch (k.first)
				{
					case 'D':
						for (std::size_t l = 0; l < k.second; ++l)
						{
							gap_histogram[pos_on_genome + l]++;
						}
					case 'M':
						pos_on_genome += k.second;

					default:
						break;
				}
			}
		}
	}

	// c. call final consensus sequence
	std::cout << "   Calling final consensus sequences" << std::endl;

	std::vector<dna_array<double, 5>> base_coverage(m_parameters.m_L, { 0.0, 0.0, 0.0, 0.0, 0.0 });
	for (auto& i : m_reads)
	{
		if (i.m_invalid == false)
		{
			int32_t pos_on_genome = i.get_opt_aln().m_POS;
			int32_t pos_in_read = i.get_opt_aln().m_left_clip_length;

			for (const auto& k : i.get_opt_aln().m_CIGAR)
			{
				switch (k.first)
				{
					case 'M':
						for (std::size_t l = 0; l < k.second; ++l)
						{
							const char base = i.m_fastq_record.m_seq[pos_in_read + l];

							if (base != 'N')
							{
								base_coverage[pos_on_genome + l][base]++;
							}
						}
						pos_in_read += k.second;
					case 'D':
						pos_on_genome += k.second;
						break;

					case 'I':
						pos_in_read += k.second;

					default:
						break;
				}
			}
		}
	}
	m_parameters.init_emission_table(std::move(base_coverage), min_freq, error_rate, ambig_bases_unequal_weight, true);

	// d. calculate edit distance and (optionally) rewrite CIGAR
	std::cout << "   Calculating edit distances";
	// replace 'M' in CIGAR if differentiate_match_state == true
	// Match: '=' Mismatch: 'X'
	if (differentiate_match_state)
	{
		std::cout << " (CIGAR 'M' will be replaced by '='/'X')";
	}
	std::cout << std::endl;

	// calculate edit/Hamming distance
	auto calculate_edit_distance = [](const reference_genome<T>& parameters, const boost::string_ref& seq, const bool rewrite_cigar, minimal_alignment& aln) -> std::tuple<uint32_t, std::string> {
		uint32_t edit_distance = 0;
		std::ostringstream MD_tag;

		uint32_t pos_in_read = aln.m_left_clip_length;
		uint32_t pos_on_genome = aln.m_POS;

		CIGAR_vec new_CIGAR_vec;

		char old_state, cur_state;
		uint32_t cur_length;

		uint32_t MD_length = 0;

		for (const auto& op : aln.m_CIGAR)
		{
			switch (op.first)
			{
				case 'M':
					if (rewrite_cigar)
					{
						old_state = (parameters.m_table_of_included_bases[pos_on_genome][seq[pos_in_read]] ? '=' : 'X');
						cur_length = 0;
					}

					for (uint32_t j = pos_in_read, i = pos_on_genome; j < pos_in_read + op.second; ++i, ++j)
					{
						cur_state = (parameters.m_table_of_included_bases[i][seq[j]] ? '=' : 'X');
						edit_distance += (cur_state == 'X');

						if (cur_state == 'X')
						{
							MD_tag << MD_length << parameters.m_ambig_ref[i];
							MD_length = 0;
						}
						else
						{
							++MD_length;
						}

						if (rewrite_cigar)
						{
							if (cur_state != old_state)
							{
								// switch states
								new_CIGAR_vec.emplace_back(old_state, cur_length);
								old_state = cur_state;
								cur_length = 0;
							}

							++cur_length;
						}
					}

					if (rewrite_cigar)
					{
						new_CIGAR_vec.emplace_back(cur_state, cur_length);
					}

					pos_on_genome += op.second;
					pos_in_read += op.second;
					break;

				case 'I':
					if (rewrite_cigar)
					{
						new_CIGAR_vec.push_back(op);
					}

					edit_distance += op.second;
					pos_in_read += op.second;
					break;

				case 'D':
					if (rewrite_cigar)
					{
						new_CIGAR_vec.push_back(op);
					}

					MD_tag << MD_length << '^' << parameters.m_ambig_ref.substr(pos_on_genome, op.second);
					MD_length = 0;

					edit_distance += op.second;
					pos_on_genome += op.second;
					break;
			}
		}

		if (rewrite_cigar)
		{
			aln.m_CIGAR.swap(new_CIGAR_vec);
		}
		MD_tag << MD_length;

		return std::tuple<uint32_t, std::string>(edit_distance, MD_tag.str());
	};

	for (auto& i : m_reads)
	{
		std::tie(i.NM, i.MD_tag) = calculate_edit_distance(m_parameters, i.m_fastq_record.m_seq, differentiate_match_state, i.get_opt_aln());
	}

	// e. sorting reads
	std::cout << "   Sorting reads" << std::endl;
	auto comp_read_entry_by_queryname = [](
		const read_entry& lhs,
		const read_entry& rhs) noexcept
	{
		return std::tie(lhs.m_fastq_record.m_id, lhs.get_opt_aln().m_POS, lhs.m_reverse_compl) < std::tie(rhs.m_fastq_record.m_id, rhs.get_opt_aln().m_POS, rhs.m_reverse_compl);
	};
	std::sort(m_reads.begin(), m_reads.end(), comp_read_entry_by_queryname);

	// f. check pairs
	std::cout << "   Setting SAM properties" << std::endl;
	m_good_reads.clear();
	m_bad_reads.clear();
	flag_reads();

	std::cout << std::endl
			  << "\tValid:                   " << m_good_reads.size() << " (" << std::fixed << std::setprecision(1) << static_cast<double>(m_good_reads.size()) * 100 / m_reads.size() << "%)" << std::endl
			  << "\tInvalid:                 " << m_bad_reads.size() << std::endl
			  << std::endl;
}

template <typename T>
void single_end_aligner<T>::flag_reads() noexcept
{
	for (auto& i : m_reads)
	{
		if (i.m_invalid == true)
		{
			m_bad_reads.push_back(&i);
		}
		else
		{
			m_good_reads.push_back(&i);
		}
	}
}

template <typename T>
void paired_end_aligner<T>::flag_reads() noexcept
{
	auto it = m_reads.begin();
	auto it_end = m_reads.end() - 1;

	std::cout << "   Performing first pass" << std::endl;
	boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> acc;

	while (it < it_end)
	{
		read_entry& first_read = *it;
		read_entry& second_read = *(it + 1);

		if (first_read.m_fastq_record.m_id == second_read.m_fastq_record.m_id)
		{
			// have paired reads
			first_read.m_paired_read = true;
			second_read.m_paired_read = true;

			first_read.m_mate_reverse_compl = second_read.m_reverse_compl;
			second_read.m_mate_reverse_compl = first_read.m_reverse_compl;

			const int32_t TLEN = std::max<decltype(TLEN)>(first_read.get_opt_aln().m_POS + first_read.get_opt_aln().m_segment_length, second_read.get_opt_aln().m_POS + second_read.get_opt_aln().m_segment_length)
				- std::min<decltype(TLEN)>(first_read.get_opt_aln().m_POS, second_read.get_opt_aln().m_POS);

			first_read.TLEN = TLEN;
			second_read.TLEN = -TLEN;

			first_read.RNEXT = '=';
			second_read.RNEXT = '=';

			first_read.PNEXT = second_read.get_opt_aln().m_POS;
			second_read.PNEXT = first_read.get_opt_aln().m_POS;

			if (first_read.m_fastq_record.m_second_in_pair + second_read.m_fastq_record.m_second_in_pair != 1)
			{
				// can only have one R1 and one R2
				first_read.m_invalid = true;
				second_read.m_invalid = true;
			}

			if ((first_read.m_reverse_compl != false) || (second_read.m_reverse_compl != true))
			{
				// Need to have
				// ---->
				//    <----
				// configuration of reads
				first_read.m_invalid = true;
				second_read.m_invalid = true;
			}

			if (first_read.m_invalid + second_read.m_invalid)
			{
				// at least one read has an issue
				first_read.m_invalid = true;
				second_read.m_invalid = true;
			}
			else
			{
				// both reads passed all checks
				acc(TLEN);
			}

			++it;
		}
		else
		{
			// first_read does not have a mate
			first_read.m_invalid = true;
		}
		++it;
	}

	std::cout << "   Performing second pass" << std::endl;

	// second pass, write output
	const int32_t MAX_TLEN = boost::accumulators::mean(acc) + 3 * std::sqrt(boost::accumulators::variance(acc));
	std::cout << std::fixed
			  << "     Mean[TLEN]: " << static_cast<int32_t>(boost::accumulators::mean(acc)) << std::endl
			  << "       sd[TLEN]: " << static_cast<int32_t>(std::sqrt(boost::accumulators::variance(acc))) << std::endl
			  << "      Max[TLEN] allowed (mean + 3*sd): " << MAX_TLEN << std::endl;

	for (auto& i : m_reads)
	{
		if ((i.m_invalid == true) || (std::abs(i.TLEN) > MAX_TLEN))
		{
			m_bad_reads.push_back(&i);
		}
		else
		{
			i.m_mapped_in_proper_pair = true;
			m_good_reads.push_back(&i);
		}
	}
}

// 7. write alignment to output
void write_header(
	std::ostream& output,
	const std::string& genome_name,
	const int32_t profile_length,
	const int32_t argc,
	const char** argv) noexcept
{
	output
		<< "@HD\tVN:1.5\tSO:queryname\n"
		<< "@SQ\tSN:" << genome_name << "\tLN:" << profile_length << '\n'
		<< "@PG\tID:" PACKAGE_NAME "\tPN:" PACKAGE_NAME "\tVN:" PACKAGE_VERSION "\tCL:" << argv[0];

	for (int i = 1; i < argc; ++i)
	{
		output << ' ' << argv[i];
	}
	output << '\n';
}

template <typename T>
void single_end_aligner<T>::write_alignment_to_file(
	const clip_mode clip,
	const std::string& consensus_name,
	const std::string& data_root,
	const std::string& output_file_name,
	const std::string& rejects_file_name) noexcept
{
	std::cout << ++m_phase << ") Writing alignment (" << (read_entry::m_clip == clip_mode::soft ? "soft" : (clip == clip_mode::hard ? "hard" : "HARD")) << " clipping)" << std::endl;

	read_entry::m_clip = clip;
	read_entry::RNAME = consensus_name;

	/////////////////////
	// write SAM files //
	/////////////////////

	// write file with good alignments
	std::ofstream output_file(output_file_name);
	write_header(output_file, consensus_name, m_parameters.m_L, m_argc, m_argv);
	for (const auto& i : m_good_reads)
	{
		output_file << *i << '\n';
	}
	output_file.close();

	// write file with failed alignments
	output_file.open(rejects_file_name, (output_file_name == rejects_file_name) ? std::ofstream::app : std::ofstream::trunc);
	if (output_file_name != rejects_file_name)
	{
		write_header(output_file, consensus_name, m_parameters.m_L, m_argc, m_argv);
	}
	for (const auto& i : m_bad_reads)
	{
		output_file << *i << '\n';
	}
	output_file.close();

	///////////////////////////////////
	// write new consensus sequences //
	///////////////////////////////////

	std::cout << "   Writing new reference sequences" << std::endl;

	// majority consensus
	output_file.open(data_root + "ref_majority.fasta");
	output_file << '>' << consensus_name << '\n'
				<< m_parameters.m_majority_ref << '\n';
	output_file.close();

	// ambiguous consensus
	output_file.open(data_root + "ref_ambig.fasta");
	output_file << '>' << consensus_name << '\n'
				<< m_parameters.m_ambig_ref << '\n';
	output_file.close();

	///////////////////
	// serialise HMM //
	///////////////////

	std::string output_hmm_filename = consensus_name;
	boost::algorithm::to_lower(output_hmm_filename);
	output_hmm_filename.append(".hmm");

	output_file.open(data_root + output_hmm_filename);
	m_parameters.save_to_file(output_file);
	output_file.close();
}

} // unnamed namespace

#endif /* NGSHMMALIGN_ALIGNER_IMPL_HPP */
