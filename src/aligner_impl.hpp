#ifndef ALIGNER_IMPL_HPP
#define ALIGNER_IMPL_HPP

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

#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <ratio>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>

#include <boost/algorithm/string.hpp>
#include <boost/progress.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <omp.h>

#include "aligner.hpp"
#include "fastq.hpp"

namespace
{

fastq_entry::fastq_entry(
	const char* id_ptr,
	std::size_t id_len,
	const char* seq_ptr,
	std::size_t seq_len,
	const char* qual_ptr,
	std::size_t qual_len)
	: m_id(id_ptr, id_len),
	  m_seq(seq_ptr, seq_len),
	  m_qual(qual_ptr, qual_len) {}

read_entry::read_entry(
	const char* id_ptr,
	std::size_t id_len,
	const char* seq_ptr,
	std::size_t seq_len,
	const char* qual_ptr,
	std::size_t qual_len)
	: m_fastq_record(
		  id_ptr,
		  id_len,
		  seq_ptr,
		  seq_len,
		  qual_ptr,
		  qual_len)
{
}

inline bool comp_read_entry_by_queryname(
	const read_entry& lhs,
	const read_entry& rhs) noexcept
{
	return (lhs.m_fastq_record.m_id < rhs.m_fastq_record.m_id);
}

template <typename T, typename V>
inline void write_seq_qual(
	std::ostream& output,
	T iter_start,
	T iter_end,
	V visitor) noexcept
{
	for (T i = iter_start; i != iter_end; ++i)
	{
		output << visitor(*i);
	}
}

void read_entry::print_sequence(std::ostream& output, bool clip_bases) const noexcept
{
	const uint32_t left_clip = (clip_bases ? m_sam_record.m_left_clip_length : 0);
	const uint32_t right_clip = (clip_bases ? m_sam_record.m_right_clip_length : 0);

	if (m_sam_record.m_forward)
	{
		// forward aligned
		auto it_start = m_fastq_record.m_seq.cbegin() + left_clip;
		auto it_end = m_fastq_record.m_seq.cend() - right_clip;
		write_seq_qual(output, it_start, it_end, identity_char);
	}
	else
	{
		// reverse complement aligned
		auto it_start = m_fastq_record.m_seq.crbegin() + left_clip;
		auto it_end = m_fastq_record.m_seq.crend() - right_clip;
		write_seq_qual(output, it_start, it_end, rev_comp_char);
	}
}

void read_entry::print_qual(std::ostream& output, bool clip_bases) const noexcept
{
	const uint32_t left_clip = (clip_bases ? m_sam_record.m_left_clip_length : 0);
	const uint32_t right_clip = (clip_bases ? m_sam_record.m_right_clip_length : 0);

	if (m_sam_record.m_forward)
	{
		// forward aligned
		auto it_start = m_fastq_record.m_qual.cbegin() + left_clip;
		auto it_end = m_fastq_record.m_qual.cend() - right_clip;
		write_seq_qual(output, it_start, it_end, identity_char);
	}
	else
	{
		// reverse complement aligned
		auto it_start = m_fastq_record.m_qual.crbegin() + left_clip;
		auto it_end = m_fastq_record.m_qual.crend() - right_clip;
		write_seq_qual(output, it_start, it_end, identity_char);
	}
}

// output
std::ostream& operator<<(
	std::ostream& output,
	const read_entry& read) noexcept
{
	auto write_clip_CIGAR = [](std::ostream& output, std::size_t length, clip_mode clip)
	{
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
	output
		<< read.m_fastq_record.m_id << '\t'
		<< read.m_sam_record.FLAG << '\t'
		<< (read.m_sam_record.RNAME ? read.m_sam_record.RNAME : "*") << '\t'
		<< read.m_sam_record.POS + 1 << '\t'
		<< read.m_sam_record.MAPQ << '\t';

	///////////
	// CIGAR //
	///////////
	write_clip_CIGAR(output, read.m_sam_record.m_left_clip_length, read.m_sam_record.m_clip);
	for (const auto& i : read.m_sam_record.CIGAR)
	{
		output << i.second << i.first;
	}
	write_clip_CIGAR(output, read.m_sam_record.m_right_clip_length, read.m_sam_record.m_clip);

	//////////////////////////
	// RNEXT + PNEXT + TLEN //
	//////////////////////////
	output
		<< '\t' << read.m_sam_record.RNEXT
		<< '\t' << read.m_sam_record.PNEXT + 1
		<< '\t' << read.m_sam_record.TLEN << '\t';

	////////////////
	// SEQ + QUAL //
	////////////////
	bool clip_bases = ((read.m_sam_record.m_clip == clip_mode::hard) || (read.m_sam_record.m_clip == clip_mode::HARD));
	read.print_sequence(output, clip_bases);
	output << '\t';
	read.print_qual(output, clip_bases);

	/////////////////////////////////////
	// Alignment Score + Edit Distance //
	/////////////////////////////////////
	return output
		<< '\t' << "NM:i:" << read.m_sam_record.NM
		<< '\t' << "MD:Z:" << read.m_sam_record.MD
		<< '\t' << "AS:i:" << read.m_sam_record.SCORE << '\n';
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

// 1. ctor
template <typename T>
std::unique_ptr<single_end_aligner<T>> single_end_aligner<T>::create_aligner_instance(
	const std::vector<std::string>& input_files,
	const int32_t min_mapped_length,
	const int argc,
	const char** argv,
	const bool write_unpaired) noexcept
{
	switch (input_files.size())
	{
		case 0:
			std::cerr << "ERROR: You have provided no input files. ngshmmalign takes either 1 (single-end) or 2 (paired-end) input file(s).\n";
			exit(EXIT_FAILURE);
			break;

		case 1:
			return std::unique_ptr<single_end_aligner<T>>(new single_end_aligner<T>(min_mapped_length, argc, argv));
			break;

		case 2:
			return std::unique_ptr<single_end_aligner<T>>(new paired_end_aligner<T>(min_mapped_length, argc, argv, write_unpaired));
			break;

		default:
			std::cerr << "ERROR: You have provided too many input files. ngshmmalign takes either 1 (single-end) or 2 (paired-end) input file(s).\n";
			exit(EXIT_FAILURE);
			break;
	}
}

template <typename T>
single_end_aligner<T>::single_end_aligner(
	const int32_t min_mapped_length_,
	const int argc_,
	const char** argv_) noexcept
	: m_argc(argc_),
	  m_argv(argv_),
	  m_min_mapped_length(min_mapped_length_)
{
}

template <typename T>
paired_end_aligner<T>::paired_end_aligner(
	const int32_t min_mapped_length_,
	const int argc_,
	const char** argv_,
	const bool write_unpaired) noexcept
	: single_end_aligner<T>(min_mapped_length_, argc_, argv_),
	  m_write_unpaired(write_unpaired)
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
	if (input_files.size() != 1)
	{
		std::terminate();
	}

	m_read_file_name = input_files.front();
	std::cout << ++m_phase << ") Loading reads into memory from " << m_read_file_name << '\n';
	m_reads_mmap.open(m_read_file_name);
	m_reads = fastq_read<read_entry>(m_reads_mmap);
}

template <typename T>
void paired_end_aligner<T>::load_reads_impl(
	const std::vector<std::string>& input_files) noexcept
{
	if (input_files.size() != 2)
	{
		std::terminate();
	}

	single_end_aligner<T>::load_reads_impl(std::vector<std::string>{ input_files.front() });

	m_read_file_name2 = input_files.back();
	std::cout << "   Loading reads into memory from " << m_read_file_name2 << '\n';
	m_reads2_mmap.open(m_read_file_name2);
	m_reads2 = fastq_read<read_entry>(m_reads2_mmap);
}

// 3. load parameters
template <typename T>
void single_end_aligner<T>::load_parameters(
	const std::string& input_file,
	background_rates& error_rates,
	const bool ambig_bases_unequal_weight) noexcept
{
	std::cout << ++m_phase << ") Loading profile HMM parameters\n";

	uint32_t L = get_length_profile();
	std::cout << "   Using Illumina length profile " << L << "\n   Parameters for pHMM:\n";

	if (error_rates.end_prob == MAGIC_NUMBER)
	{
		error_rates.end_prob = 1.0 / L;
	}
	if (error_rates.right_clip_open == MAGIC_NUMBER)
	{
		error_rates.right_clip_open = error_rates.left_clip_open / L;
	}
	if (m_min_mapped_length == std::numeric_limits<decltype(m_min_mapped_length)>::max())
	{
		m_min_mapped_length = (L * 4) / 5;
	}
	std::cout << '\n' << error_rates << '\t' << "Minimum mapped length:   " << m_min_mapped_length << "\n\n";

	std::string ref_ext(boost::filesystem::extension(input_file));
	if ((ref_ext == ".fasta") || (ref_ext == ".fas") || (ref_ext == ".fna"))
	{
		std::cout << "   Input = FASTA\n";
		m_parameters.set_parameters(
			input_file,
			error_rates,
			L,
			ambig_bases_unequal_weight,
			typename reference_genome<T>::msa_input{});
	}
	else if (ref_ext == ".hmm")
	{
		std::cout << "   Input = " PACKAGE " serialised format\n";
		m_parameters.set_parameters(
			input_file,
			error_rates,
			L,
			ambig_bases_unequal_weight,
			typename reference_genome<T>::serialized_input{});
	}
	else
	{
		std::cerr << "ERROR: file ending" << ref_ext << " not recognized\n";
		exit(EXIT_FAILURE);
	}

	if (m_parameters.m_reference_genome_name.empty())
	{
		// have a multiple sequence alignment
		serialize = true;
	}
	else
	{
		// only a single reference sequence given
		serialize = false;
	}
}

// 4. sort reads
template <typename T>
void single_end_aligner<T>::sort_reads() noexcept
{
	std::cout << ++m_phase << ") Sorting reads from " << m_read_file_name << '\n';
	sort_reads_impl();
}

template <typename T>
void single_end_aligner<T>::sort_reads_impl() noexcept
{
	std::sort(m_reads.begin(), m_reads.end(), comp_read_entry_by_queryname);
}

template <typename T>
void paired_end_aligner<T>::sort_reads_impl() noexcept
{
	single_end_aligner<T>::sort_reads_impl();

	std::cout << "   Sorting reads from " << m_read_file_name2 << '\n';
	std::sort(m_reads2.begin(), m_reads2.end(), comp_read_entry_by_queryname);
}

// 5. perform parameter estimation
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

	partioned_genome() = delete;
	partioned_genome(const partioned_genome&) = delete;
	partioned_genome(partioned_genome&&) = delete;
	partioned_genome& operator=(const partioned_genome& other) = delete;
	partioned_genome& operator=(partioned_genome&& other) = delete;

	partioned_genome(const std::string& tmpdir, const int32_t L, const int32_t read_profile) noexcept
		: m_offset(read_profile / 6),
		  m_num_bins(L / m_offset)
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
		//uint64_t num_unique_assignments = 0;
		//uint64_t num_multiple_assignments = 0;
		//uint64_t num_no_assignments = 0;
		const int32_t max_index = m_regions.size() - 1;

		for (const auto& read : reads)
		{
			const auto start = read.m_sam_record.POS - read.m_sam_record.m_left_clip_length;
			const auto end = read.m_sam_record.POS + read.m_sam_record.m_segment_length + read.m_sam_record.m_right_clip_length;

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

				//if (candidate_regions.size() == 1)
				//{
				//	++num_unique_assignments;
				//}
				//else
				//{
				//	++num_multiple_assignments;
				//}
			}
			//else
			//{
			//	++num_no_assignments;
			//}
		}

		//std::cout
		//	<< "Unique:   " << num_unique_assignments << '\n'
		//	<< "Multiple: " << num_multiple_assignments << '\n'
		//	<< "None:     " << num_no_assignments << '\n';
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
				i.m_reads[j]->print_sequence(output, false);
				output << '\n';
			}
			output << ">REF\n" << reference.substr(i.m_start, i.m_end - i.m_start) << '\n';
			output.close();
		}
	}

	void run_mafft_on_files(const std::string& mafft, const int32_t num_threads) const noexcept
	{
		std::cout << "Running MAFFT on region\n";
		for (const auto& i : m_regions)
		{
			std::ostringstream mafft_cmd_line;
			mafft_cmd_line
				<< mafft
				<< " --thread " << num_threads
				<< " --localpair --maxiterate 1000 --quiet "
				<< i.m_files_dir << "/reads.fasta > "
				<< i.m_files_dir << "/reads_aln.fasta";

			std::cout << std::right << std::setw(8) << i.m_start << " -> " << std::setw(6) << i.m_end << '\n';
			int result = std::system(mafft_cmd_line.str().c_str());

			if (result != 0)
			{
				std::cerr << "ERROR: MAFFT exited with error code " << result << '\n';
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

		std::cout << "Processing MAFFT aligned files\n" << std::fixed;
		for (const auto& i : m_regions)
		{
			std::cout << std::right << std::setw(8) << i.m_start << " -> " << std::setw(6) << i.m_end << '\n';

			// 0. load aligned MSA
			std::vector<reference_haplotype> msa_reads(fasta_read<reference_haplotype>(i.m_files_dir + "/reads_aln.fasta"));
			if (msa_reads.empty())
			{
				// no reads in alignment
				msa_reads = fasta_read<reference_haplotype>(i.m_files_dir + "/reads.fasta");
			}
			const reference_haplotype ref_seq(msa_reads.back());

			if (msa_reads.size() > 1)
			{
				msa_reads.pop_back();
			}

			if (ref_seq.name != "REF")
			{
				std::cerr << "ERROR: Last sequence in " << i.m_files_dir << "/reads_aln.fasta is not the reference!\n";
				exit(EXIT_FAILURE);
			}

			// 1. collect statistics
			const int32_t total_length_of_MSA = ref_seq.sequence.length();
			std::vector<uint32_t> coverage(total_length_of_MSA, 0);
			std::vector<uint32_t> nongap_coverage(total_length_of_MSA, 0);
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
							++nongap_coverage[j];
							++base_dist[j][base];
						}
					}
				}
			}

			// 2. determine whether a locus is of high coverage
			std::vector<uint8_t> is_locus_cov_high;
			is_locus_cov_high.reserve(total_length_of_MSA);
			for (int32_t j = 0; j < total_length_of_MSA; ++j)
			{
				is_locus_cov_high.push_back((static_cast<double>(nongap_coverage[j]) / coverage[j]) > 0.5);
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
				// a general problem of our estimation scheme is that is rather sensitive
				// to artefacts in low coverage regimes. In order to decide whether
				// a locus in the alignment is 'valid', we therefore require at least
				// 3 reads supporting the locus AND a minimum fraction of 'min_base_cov',
				// and if we have less than 3 reads, all of them must support the position.
				if ((static_cast<double>(nongap_coverage[j]) / coverage[j]) >= std::min<double>(std::max<double>(2.5 / coverage[j], min_base_cov), 1.0))
				{
					valid_loci.push_back(j);
				}
			}

			// 5. determine parameters
			auto remove_low_frequency_bases = [min_base_cutoff](dna_array<double, 5>& allel_freq) -> void
			{
				const double coverage = std::accumulate(&allel_freq[static_cast<std::size_t>(0)], &allel_freq[static_cast<std::size_t>(4)], 0.0);
				for (char k : { 'A', 'C', 'G', 'T' })
				{
					if (allel_freq[k] / coverage < min_base_cutoff)
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
				if ((num_M_D / num_M < 0.05) || (num_M < 10))
				{
					num_M_D = 0;
				}

				num_D += (!num_D);
				if ((num_D_D / num_D < 0.05) || (num_D < 10))
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
	const bool ambig_bases_unequal_weight) noexcept
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
	std::cout << ++m_phase << ") Performing first rough alignment\n";
	perform_alignment_impl(clip_mode::soft, seed, false, verbose, false);

	// 2. chop up genome and extract reads
	partioned_genome separated_reads(tmpdir_root, m_parameters.m_L, m_parameters.read_length_profile);
	estimate_parameters_impl(separated_reads, rng);

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

	m_parameters.m_reference_genome_name.clear();
	serialize = true;
	m_parameters.set_parameters(std::move(E_p), std::move(M_D_p), std::move(D_D_p), error_rates, ambig_bases_unequal_weight);

	// 7. (optional) cleanup temporary MAFFT files
	if (keep_mafft_files == false)
	{
		boost::filesystem::remove_all(tmpdir_root);
	}
}

template <typename T>
void single_end_aligner<T>::estimate_parameters_impl(partioned_genome& separated_reads, std::default_random_engine& generator) noexcept
{
	separated_reads.add_reads(m_reads, generator);
}

template <typename T>
void paired_end_aligner<T>::estimate_parameters_impl(partioned_genome& separated_reads, std::default_random_engine& generator) noexcept
{
	single_end_aligner<T>::estimate_parameters_impl(separated_reads, generator);
	separated_reads.add_reads(m_reads2, generator);
}

// 6. perform alignment
template <typename T>
std::size_t single_end_aligner<T>::number_of_reads() const noexcept
{
	return m_reads.size();
}

template <typename T>
std::size_t paired_end_aligner<T>::number_of_reads() const noexcept
{
	return single_end_aligner<T>::number_of_reads() + m_reads2.size();
}

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
void alignment_impl(
	std::vector<read_entry>& reads,
	const reference_genome<T>& parameters,
	const clip_mode clip,
	const uint64_t seed,
	const bool exhaustive,
	const bool verbose,
	const bool differentiate_match_state) noexcept
{
	std::cerr << "T has to be an integral type!\n";
	exit(EXIT_FAILURE);
}

void display_progress(
	const unsigned long& progress,
	const unsigned long total) noexcept
{
	boost::progress_display show_progress(total);
	unsigned long now = 0, previous = 0;

	do
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(250));
		show_progress += (now - previous);

		previous = now;
		now = progress;
	} while (now < total);
}

template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
void alignment_impl(
	std::vector<read_entry>& reads,
	const reference_genome<T>& parameters,
	const clip_mode clip,
	const uint64_t seed,
	const bool exhaustive,
	const bool verbose,
	const bool differentiate_match_state) noexcept
{
	const auto end = reads.cend();
	unsigned long progress = 0;

	// minimum number of matched k-mers of a read to be considered the correct strand
	constexpr uint32_t min_hash_pos_better = 20;

	// maximum number of matched k-mers of a read to be considered the correct strand
	constexpr uint32_t max_hash_pos_worse = 2;

#ifndef NDEBUG
	struct alignment_stats_struct
	{
		uint64_t num_exhaustive_alignments = 0;
		uint64_t num_indexed_alignments = 0;
		uint64_t num_sacrificial_alignments = 0;
	} alignment_stats;
	std::cout << '\n';
#endif

	std::function<void()> display_func;
	if (verbose)
	{
#ifdef NDEBUG
		display_func = std::bind(display_progress, std::cref(progress), reads.size());
	}
	else
	{
#endif
		display_func = []()
		{
		};
	}
	std::thread display_thread(display_func);

	std::random_device rd;
#pragma omp parallel shared(reads, progress)
	{
		std::default_random_engine rng;
#pragma omp critical
		{
			uint64_t thread_seed = (seed ? seed : std::uniform_int_distribution<uint64_t>()(rd));
#ifndef NDEBUG
			int tid = omp_get_thread_num();
			std::cout << "Seeding Thread #" << tid << " with seed " << thread_seed << '\n';
#endif
			rng.seed(thread_seed);
		}

		sam_entry forward_sam, reverse_sam;
		std::string rev_seq;

#pragma omp for schedule(runtime)
		for (auto it = reads.begin(); it < end; ++it)
		{
			read_entry& read_entry_ref = *it;
			const boost::string_ref seq = read_entry_ref.m_fastq_record.m_seq;

			if (verbose)
			{
#pragma omp atomic
				++progress;
			}

			// 1. first use heuristics to determine strand
			bool align_forward = false;
			forward_sam.SCORE = std::numeric_limits<decltype(forward_sam.SCORE)>::min();
			auto forward_stats = parameters.find_pos(seq);

			bool align_reverse = false;
			reverse_sam.SCORE = std::numeric_limits<decltype(reverse_sam.SCORE)>::min();

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
				&parameters,
				&rng](sam_entry& sam_alignment, const boost::string_ref& query, bool exhaustive, const typename reference_genome<T>::index_stat& heuristics, bool differentiate_match_state)
			{
				int32_t POS;
				int32_t end_POS;

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

// // test for detection of failed heuristic
// POS = heuristics.POS + 100;
// POS = std::max<int32_t>(0, POS);
//
// end_POS = heuristics.POS + query.length() - 100;
// end_POS = std::min<int32_t>(end_POS, parameters.m_L);
#ifndef NDEBUG
					++alignment_stats.num_indexed_alignments;
#endif
				}

				sam_alignment = hmmalign<T>::viterbi(parameters, query, rng, POS, end_POS, differentiate_match_state);

#ifndef NDEBUG
				{
					std::cout
						<< "POS:" << std::right
						<< std::setw(6) << POS
						<< std::setw(6) << "End:"
						<< std::setw(6) << end_POS
						<< std::setw(14) << "True POS:"
						<< std::setw(6) << sam_alignment.POS
						<< std::setw(14) << "True End:"
						<< std::setw(6) << sam_alignment.POS + sam_alignment.m_segment_length << '\n';
				}
#endif

				if ((((sam_alignment.POS - sam_alignment.m_left_clip_length < POS) && (POS != 0)) || ((sam_alignment.POS + sam_alignment.m_segment_length + sam_alignment.m_right_clip_length > end_POS) && (end_POS != static_cast<int32_t>(parameters.m_L))))
					&& (exhaustive == false))
				{
					// window size probably too small, perform sacrificial full-length alignment
					sam_alignment = hmmalign<T>::viterbi(parameters, query, rng, 0, parameters.m_L, differentiate_match_state);
#ifndef NDEBUG
					{
						++alignment_stats.num_sacrificial_alignments;
						std::cout << "Had to perform sacrificial exhaustive alignment\n";
					}
#endif
				}
			};

			if (align_forward)
			{
				per_strand_alignment(forward_sam, seq, exhaustive_alignment, forward_stats, differentiate_match_state);
			}

			if (align_reverse)
			{
				per_strand_alignment(reverse_sam, rev_seq, exhaustive_alignment, reverse_stats, differentiate_match_state);
			}

			// 3. pick correct alignment
			if (reverse_sam.SCORE > forward_sam.SCORE)
			{
				// reverse complementary is better
				read_entry_ref.m_sam_record.m_forward = false;

				read_entry_ref.m_sam_record = std::move(reverse_sam);
				read_entry_ref.m_sam_record.FLAG |= 0x10;
				read_entry_ref.m_sam_record.m_forward = false;
			}
			else
			{
				// forward is better
				read_entry_ref.m_sam_record.m_forward = true;

				read_entry_ref.m_sam_record = std::move(forward_sam);
			}
			read_entry_ref.m_sam_record.m_clip = clip;
		}
	}

#ifndef NDEBUG
	std::cout
		<< "Number of exhaustive alignments:  " << alignment_stats.num_exhaustive_alignments << '\n'
		<< "Number of indexed alignments:     " << alignment_stats.num_indexed_alignments << '\n'
		<< "Number of sacrificial alignments: " << alignment_stats.num_sacrificial_alignments << '\n';
#endif

	display_thread.join();
}

template <typename T>
void single_end_aligner<T>::perform_alignment(
	const std::string& reference_genome_name,
	const clip_mode clip,
	const uint64_t seed,
	const bool exhaustive,
	const bool verbose,
	const bool differentiate_match_state) noexcept
{
	std::cout << ++m_phase << ") Performing alignment (" << (clip == clip_mode::soft ? "soft" : (clip == clip_mode::hard ? "hard" : "HARD")) << " clipping)\n";
	if (m_parameters.m_reference_genome_name.empty())
	{
		m_parameters.m_reference_genome_name = reference_genome_name;
	}
	if (seed)
	{
		std::cout << "   Using provided seed " << seed << " for deterministic alignment\n";
	}
	if (differentiate_match_state)
	{
		std::cout << "   CIGAR 'M' will be replaced by '='/'X'\n";
	}
	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	perform_alignment_impl(clip, seed, exhaustive, verbose, differentiate_match_state);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	const double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1E9;
	const double performance = number_of_reads() / duration;

	std::cout
		<< "   Alignment took "
		<< static_cast<uint32_t>(duration) / 3600 << "h "
		<< (static_cast<uint32_t>(duration) / 60) % 60 << "min "
		<< static_cast<uint32_t>(duration) % 60 << "s ("
		<< std::fixed << std::setprecision(1)
		<< performance / num_threads << " reads/thread/s, " << performance << " reads/s).\n";
}

template <typename T>
void single_end_aligner<T>::perform_alignment_impl(
	const clip_mode clip,
	const uint64_t seed,
	const bool exhaustive,
	const bool verbose,
	const bool differentiate_match_state) noexcept
{
	std::cout << "   Aligning reads from " << m_read_file_name << " (" << m_reads.size() << ")" << std::flush;
	alignment_impl(m_reads, m_parameters, clip, seed, exhaustive, verbose, differentiate_match_state);
	std::cout << '\n';
}

template <typename T>
void paired_end_aligner<T>::perform_alignment_impl(
	const clip_mode clip,
	const uint64_t seed,
	const bool exhaustive,
	const bool verbose,
	const bool differentiate_match_state) noexcept
{
	single_end_aligner<T>::perform_alignment_impl(clip, seed, exhaustive, verbose, differentiate_match_state);

	std::cout << "   Aligning reads from " << m_read_file_name2 << " (" << m_reads2.size() << ")" << std::flush;
	alignment_impl(m_reads2, m_parameters, clip, seed, exhaustive, verbose, differentiate_match_state);
	std::cout << '\n';
}

// 7. write alignment to output
template <typename T>
void single_end_aligner<T>::write_alignment_to_file(
	const std::string& data_root,
	const std::string& output_file_name,
	const std::string& rejects_file_name) noexcept
{
	std::cout << ++m_phase << ") Writing ";
	write_alignment_to_file_impl(output_file_name, rejects_file_name);
	std::cout << '\n';

	if (serialize)
	{
		// provided an MSA
		std::cout << "   Writing new reference sequences\n";

		std::ofstream output;
		output.open(data_root + "ref_majority.fasta");
		output << '>' << m_parameters.m_reference_genome_name << '\n' << m_parameters.m_majority_ref << '\n';
		output.close();

		output.open(data_root + "ref_ambig.fasta");
		output << '>' << m_parameters.m_reference_genome_name << '\n' << m_parameters.m_ambig_ref << '\n';
		output.close();

		std::string output_hmm_filename = m_parameters.m_reference_genome_name;
		boost::algorithm::to_lower(output_hmm_filename);
		output_hmm_filename.append(".hmm");

		output.open(data_root + output_hmm_filename);
		m_parameters.save_to_file(output);
		output.close();
	}
}

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
void single_end_aligner<T>::write_alignment_to_file_impl(
	const std::string& output_file_name,
	const std::string& rejects_file_name) noexcept
{
	std::cout << "single-end alignment\n";
	std::size_t proper_reads = 0;
	std::size_t too_short_reads = 0;

	std::ofstream output_file(output_file_name);
	std::ofstream rejects_file(rejects_file_name);

	write_header(output_file, m_parameters.m_reference_genome_name, m_parameters.m_L, m_argc, m_argv);
	if (output_file_name != rejects_file_name)
	{
		write_header(rejects_file, m_parameters.m_reference_genome_name, m_parameters.m_L, m_argc, m_argv);
	}

	for (const auto& i : m_reads)
	{
		((i.m_sam_record.m_mapped_length < m_min_mapped_length) ? (++too_short_reads, rejects_file) : (++proper_reads, output_file)) << i;
	}

	std::cout
		<< "\n\tToo short:               " << too_short_reads << '\n'
		<< '\n'
		<< "\tReads                    " << proper_reads << " (" << std::fixed << std::setprecision(1) << static_cast<double>(proper_reads) * 100 / m_reads.size() << "%)\n";
}

template <typename T>
void paired_end_aligner<T>::write_alignment_to_file_impl(
	const std::string& output_file_name,
	const std::string& rejects_file_name) noexcept
{
	std::cout << "paired-end alignment\n";

	/* first pass, set and correct flags and collect statistics */
	const auto end1 = m_reads.end();
	const auto end2 = m_reads2.end();

	auto it1 = m_reads.begin();
	auto it2 = m_reads2.begin();

	int32_t TLEN;
	boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> acc;

	for (; (it1 != end1) && (it2 != end2);)
	{
		// check whether reads are properly paired
		if (comp_read_entry_by_queryname(*it1, *it2))
		{
			++it1;
		}
		else
		{
			if (!(comp_read_entry_by_queryname(*it2, *it1)))
			{
				// swap reads depending on position and other criteria
				if (it2->m_sam_record.POS < it1->m_sam_record.POS)
				{
					std::swap(*it1, *it2);
				}
				else
				{
					if (it2->m_sam_record.POS == it1->m_sam_record.POS)
					{
						if ((it2->m_sam_record.m_forward) && (it1->m_sam_record.m_forward == false))
						{
							std::swap(*it1, *it2);
						}
						else
						{
							if (it2->m_fastq_record.m_seq < it1->m_fastq_record.m_seq)
							{
								std::swap(*it1, *it2);
							}
						}
					}
				}

				it1->m_sam_record.FLAG |= 0x3;
				it2->m_sam_record.FLAG |= 0x3;

				it1->m_sam_record.FLAG |= 0x40;
				it2->m_sam_record.FLAG |= 0x80;

				if (it1->m_sam_record.m_forward == false)
				{
					it2->m_sam_record.FLAG |= 0x20;
				}
				if (it2->m_sam_record.m_forward == false)
				{
					it1->m_sam_record.FLAG |= 0x20;
				}

				it1->m_sam_record.RNEXT = '=';
				it2->m_sam_record.RNEXT = '=';

				it1->m_sam_record.PNEXT = it2->m_sam_record.POS;
				it2->m_sam_record.PNEXT = it1->m_sam_record.POS;

				TLEN = (it2->m_sam_record.POS + it2->m_sam_record.m_segment_length) - (it1->m_sam_record.POS);
				it1->m_sam_record.TLEN = TLEN;
				it2->m_sam_record.TLEN = -TLEN;

				acc(TLEN);
				++it1;
			}

			++it2;
		}
	}

	/* second pass, write output */
	const int32_t MAX_TLEN = boost::accumulators::mean(acc) + 3 * std::sqrt(boost::accumulators::variance(acc));
	std::cout << std::fixed
			  << "     Mean[TLEN]: " << static_cast<int32_t>(boost::accumulators::mean(acc))
			  << "\n       sd[TLEN]: " << static_cast<int32_t>(std::sqrt(boost::accumulators::variance(acc)))
			  << "\n      Max[TLEN] allowed (mean + 3*sd): " << MAX_TLEN << '\n';

	std::size_t unpaired_reads = 0;
	std::size_t too_short_reads = 0;
	std::size_t TLEN_too_big = 0;

	std::size_t num_for_for = 0;
	std::size_t num_rev_for = 0;
	std::size_t num_for_rev = 0;
	std::size_t num_rev_rev = 0;

	std::ofstream output_file(output_file_name);
	std::ofstream rejects_file(rejects_file_name);
	std::ofstream& unpaired_file = (m_write_unpaired ? output_file : rejects_file);

	write_header(output_file, m_parameters.m_reference_genome_name, m_parameters.m_L, m_argc, m_argv);
	if (output_file_name != rejects_file_name)
	{
		write_header(rejects_file, m_parameters.m_reference_genome_name, m_parameters.m_L, m_argc, m_argv);
	}

	it1 = m_reads.begin();
	it2 = m_reads2.begin();

	for (; (it1 != end1) && (it2 != end2);)
	{
		// check whether reads are too short
		if (it1->m_sam_record.m_mapped_length < m_min_mapped_length)
		{
			++too_short_reads;
			rejects_file << *it1;
			++it1;
			continue;
		}
		if (it2->m_sam_record.m_mapped_length < m_min_mapped_length)
		{
			++too_short_reads;
			rejects_file << *it2;
			++it2;
			continue;
		}

		// check whether reads are properly paired
		if (comp_read_entry_by_queryname(*it1, *it2))
		{
			++unpaired_reads;
			unpaired_file << *it1;
			++it1;
		}
		else
		{
			if (!(comp_read_entry_by_queryname(*it2, *it1)))
			{
				if (it1->m_sam_record.TLEN < MAX_TLEN)
				{
					switch (it1->m_sam_record.FLAG)
					{
						case 67: // + 131
							// ---> ... --->
							assert(it2->m_sam_record.FLAG == 131);
							num_for_for += 2;
							rejects_file << *it1;
							rejects_file << *it2;
							break;

						case 83: // + 163
							// <--- ... --->
							assert(it2->m_sam_record.FLAG == 163);
							num_rev_for += 2;
							rejects_file << *it1;
							rejects_file << *it2;
							break;

						case 99: // + 147
							// ---> ... <---
							assert(it2->m_sam_record.FLAG == 147);
							num_for_rev += 2;
							output_file << *it1;
							output_file << *it2;
							break;

						case 115: // + 179
							// <--- ... <---
							assert(it2->m_sam_record.FLAG == 179);
							num_rev_rev += 2;
							rejects_file << *it1;
							rejects_file << *it2;
							break;

						default:
							std::cerr << it1->m_sam_record.FLAG << " is an unrecognised flag in paired-end alignments\n";
							exit(EXIT_FAILURE);
					}
				}
				else
				{
					TLEN_too_big += 2;
					rejects_file << *it1;
					rejects_file << *it2;
				}

				++it1;
			}
			else
			{
				++unpaired_reads;
				unpaired_file << *it2;
			}

			++it2;
		}
	}

	for (; it1 != end1; ++it1)
	{
		(it1->m_sam_record.m_mapped_length < m_min_mapped_length ? (++too_short_reads, rejects_file) : (++unpaired_reads, unpaired_file)) << *it1;
	}
	for (; it2 != end2; ++it2)
	{
		(it2->m_sam_record.m_mapped_length < m_min_mapped_length ? (++too_short_reads, rejects_file) : (++unpaired_reads, unpaired_file)) << *it2;
	}

	std::cout
		<< "\n\tToo short:               " << too_short_reads << '\n'
		<< "\t" << (m_write_unpaired ? "Keeping unpaired:   " : "Discarding unpaired:") << "     " << unpaired_reads << '\n'
		<< "\tTLEN/Fragment too large: " << TLEN_too_big << '\n'
		<< "\tReads ---> ... --->:     " << num_for_for << '\n'
		<< "\tReads <--- ... --->:     " << num_rev_for << '\n'
		<< "\tReads <--- ... <---:     " << num_rev_rev << '\n'
		<< '\n'
		<< "\tReads ---> ... <---:     " << num_for_rev << " (" << std::fixed << std::setprecision(1) << static_cast<double>(num_for_rev) * 100 / (m_reads.size() + m_reads2.size()) << "%)\n";
}
}

#endif /* ALIGNER_IMPL_HPP */