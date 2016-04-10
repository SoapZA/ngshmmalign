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

#include <boost/progress.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <omp.h>

#include "aligner.hpp"
#include "fastq.hpp"

/* fastq_entry */
template <typename T1, typename T2, typename T3>
fastq_entry::fastq_entry(T1&& id_, T2&& seq_, T3&& qual_)
	: m_id(std::forward<T1>(id_)), m_seq(std::forward<T2>(seq_)), m_qual(std::forward<T3>(qual_))
{
	static_assert(std::is_same<T1, std::string>::value, "T1 has to be of type std::string!\n");
	static_assert(std::is_same<T2, std::string>::value, "T2 has to be of type std::string!\n");
	static_assert(std::is_same<T3, std::string>::value, "T3 has to be of type std::string!\n");
}

/* read_entry */
template <typename T1, typename T2, typename T3>
read_entry::read_entry(T1&& id_, T2&& seq_, T3&& qual_)
	: m_fastq_record(std::forward<T1>(id_), std::forward<T2>(seq_), std::forward<T3>(qual_))
{
	static_assert(std::is_same<T1, std::string>::value, "T1 has to be of type std::string!\n");
	static_assert(std::is_same<T2, std::string>::value, "T2 has to be of type std::string!\n");
	static_assert(std::is_same<T3, std::string>::value, "T3 has to be of type std::string!\n");
}

inline bool read_entry::operator<(const read_entry& other_read_entry) const noexcept
{
	return (m_fastq_record.m_id < other_read_entry.m_fastq_record.m_id);
}

// output
std::ostream& operator<<(std::ostream& output, const read_entry& read) noexcept
{
	return output
		<< read.m_fastq_record.m_id << '\t'
		<< read.m_sam_record.FLAG << '\t'
		<< "CONSENSUS" << '\t'
		<< read.m_sam_record.POS << '\t'
		<< read.m_sam_record.MAPQ << '\t'
		<< read.m_sam_record.CIGAR << '\t'
		<< read.m_sam_record.RNEXT << '\t'
		<< read.m_sam_record.PNEXT << '\t'
		<< read.m_sam_record.TLEN << '\t'
		<< read.m_fastq_record.m_seq << '\t'
		<< read.m_fastq_record.m_qual << '\t'
		<< "AS:i:" << read.m_sam_record.Score << '\n';
}

template <typename T>
uint32_t single_end_aligner<T>::get_length_profile() const noexcept
{
	int32_t max_length = 0;
	for (std::vector<read_entry>::size_type i = 0; i < std::min(10000ul, m_reads.size()); ++i)
	{
		max_length = std::max<int32_t>(max_length, m_reads[i].m_fastq_record.m_seq.length());
	}

	int32_t temp, diff = std::numeric_limits<int32_t>::max();
	uint32_t result;

	for (auto i : { 25, 36, 50, 75, 150, 250, 300 }) // Illumina profiles
	{
		temp = abs(max_length - i);
		if (temp < diff)
		{
			diff = temp;
			result = i;
		}
	}

	return result;
}

/* single_end_aligner */
template <typename T>
single_end_aligner<T>::single_end_aligner(const std::string& file_name_, uint32_t min_mapped_length_, int argc_, char** argv_) noexcept
	: m_argc(argc_),
	  m_argv(argv_),
	  m_min_mapped_length(min_mapped_length_),
	  m_file_name(file_name_),
	  m_reads(fastq_read<read_entry>(m_file_name))
{
	std::cout << "2) Loading reads into memory from " << m_file_name << '\n';
}

template <typename T>
paired_end_aligner<T>::paired_end_aligner(const std::string& file_name_, const std::string& file_name2_, bool write_unpaired, uint32_t min_mapped_length_, int argc_, char** argv_)
	: single_end_aligner<T>(file_name_, min_mapped_length_, argc_, argv_), m_file_name2(file_name2_), m_reads2(fastq_read<read_entry>(m_file_name2)), m_write_unpaired(write_unpaired)
{
	std::cout << "   Loading reads into memory from " << m_file_name2 << '\n';
}

// mutators
template <typename T>
template <typename V>
void single_end_aligner<T>::set_parameters(
	const std::vector<dna_array<V, 5>>& allel_freq_,
	const std::vector<V>& vec_M_to_D_p_,
	const std::vector<V>& vec_D_to_D_p_,
	background_rates<V> error_rates) noexcept
{
	uint32_t L = get_length_profile();
	std::cout << "4) Using Illumina length profile " << L << "\n   Parameters for pHMM:\n";

	if (error_rates.end_prob == MAGIC_NUMBER)
	{
		error_rates.end_prob = 1.0 / L;
	}

	if (error_rates.right_clip_open == MAGIC_NUMBER)
	{
		error_rates.right_clip_open = error_rates.left_clip_open / L;
	}

	if (m_min_mapped_length == MAGIC_NUMBER)
	{
		m_min_mapped_length = (L * 4) / 5;
	}

	std::cout << '\n' << error_rates << '\t' << "Minimum mapped length:   " << m_min_mapped_length << "\n\n";

	m_profile_hmm.set_parameters(
		allel_freq_,
		vec_M_to_D_p_,
		vec_D_to_D_p_,
		error_rates);
}

// sort reads
template <typename T>
void single_end_aligner<T>::sort_reads() noexcept
{
	std::cout << "3) Sorting reads from " << m_file_name << '\n';
	sort_reads_impl();
}

template <typename T>
void single_end_aligner<T>::sort_reads_impl() noexcept
{
	std::sort(m_reads.begin(), m_reads.end());
}

template <typename T>
void paired_end_aligner<T>::sort_reads_impl() noexcept
{
	single_end_aligner<T>::sort_reads_impl();

	std::cout << "   Sorting reads from " << m_file_name2 << '\n';
	std::sort(m_reads2.begin(), m_reads2.end());
}

// number of reads
template <typename T>
std::size_t single_end_aligner<T>::number_of_reads() const noexcept
{
	return number_of_reads_impl();
}

template <typename T>
std::size_t single_end_aligner<T>::number_of_reads_impl() const noexcept
{
	return m_reads.size();
}

template <typename T>
std::size_t paired_end_aligner<T>::number_of_reads_impl() const noexcept
{
	return single_end_aligner<T>::number_of_reads_impl() + m_reads2.size();
}

// perform alignment
inline std::string rev_comp(const std::string& str) noexcept
{
	std::string result;
	result.resize(str.length());

	std::transform(str.rbegin(), str.rend(), result.begin(),
		[](char base) -> char
		{
			switch (base)
			{
				case 'A':
					return 'T';
				case 'C':
					return 'G';
				case 'G':
					return 'C';
				case 'T':
					return 'A';
				case 'N':
					return 'N';
				default:
					throw std::invalid_argument(std::string(1, base) + std::string(" is an unknown base! Exiting...\n"));
			}
		});
	return result;
}

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
void alignment_impl(std::vector<read_entry>& reads, const hmmalign<T>& profile_hmm, char clip_char, bool verbose) noexcept
{
	std::cerr << "T has to be an integral type!\n";
	exit(EXIT_FAILURE);
}

uint64_t get_random_seed()
{
	uint64_t random_seed;
	std::ifstream file("/dev/urandom", std::ios::binary);
	if (file.is_open())
	{
		char* memblock;
		uint64_t size = sizeof(uint64_t);
		memblock = new char[size];
		file.read(memblock, size);
		file.close();
		random_seed = *reinterpret_cast<uint64_t*>(memblock);
		delete[] memblock;
	}
	else
	{
		std::cerr << "Could not open /dev/urandom for generating a seed!\n";
		exit(EXIT_FAILURE);
	}

	return random_seed;
}

void display_progress(const unsigned long& progress, const unsigned long total)
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
void alignment_impl(std::vector<read_entry>& reads, const hmmalign<T>& profile_hmm, char clip_char, bool verbose) noexcept
{
	const auto end = reads.end();
	unsigned long progress = 0;

	std::function<void()> display_func;
	if (verbose)
	{
		display_func = std::bind(display_progress, std::cref(progress), reads.size());
	}
	else
	{
		display_func = [](){};
	}
	std::thread display_thread(display_func);

#pragma omp parallel shared(reads, progress)
	{
		std::default_random_engine rng(get_random_seed());

		sam_entry forward_sam, reverse_sam;
		std::string rev_seq;

#pragma omp for schedule(dynamic, 200)
		for (auto it = reads.begin(); it < end; ++it)
		{
			read_entry& read_entry_ref = *it;
			std::string& seq = read_entry_ref.m_fastq_record.m_seq;

			if (verbose)
			{
#pragma omp atomic
				++progress;
			}

			// 1. try forward
			sam_entry forward_sam = profile_hmm.viterbi(seq, rng, clip_char);

			// 2. try reverse
			rev_seq = rev_comp(seq);
			sam_entry reverse_sam = profile_hmm.viterbi(rev_seq, rng, clip_char);

			// 3. pick correct alignment
			if (reverse_sam.Score > forward_sam.Score)
			{
				// reverse complementary is better
				read_entry_ref.m_fastq_record.m_seq = std::move(rev_seq);
				std::reverse(read_entry_ref.m_fastq_record.m_qual.begin(), read_entry_ref.m_fastq_record.m_qual.end());

				read_entry_ref.m_sam_record = std::move(reverse_sam);
				read_entry_ref.m_sam_record.FLAG |= 0x10;
				read_entry_ref.m_sam_record.forward = false;
			}
			else
			{
				// forward is better
				read_entry_ref.m_sam_record = std::move(forward_sam);
			}

			// 4. Clip bases (if necessary)
			if (clip_char == 'H')
			{
				std::string& seq = read_entry_ref.m_fastq_record.m_seq;
				seq = seq.substr(read_entry_ref.m_sam_record.Left_clip_length, seq.length() - read_entry_ref.m_sam_record.Left_clip_length - read_entry_ref.m_sam_record.Right_clip_length);

				std::string& qual = read_entry_ref.m_fastq_record.m_qual;
				qual = qual.substr(read_entry_ref.m_sam_record.Left_clip_length, qual.length() - read_entry_ref.m_sam_record.Left_clip_length - read_entry_ref.m_sam_record.Right_clip_length);
			}
		}
	}

	display_thread.join();
}

template <typename T>
void single_end_aligner<T>::perform_alignment(char clip_char, bool verbose) noexcept
{
	std::cout << "5) Performing alignment (" << (clip_char == 'H' ? "Hard" : "Soft") << " clipping) with -t = " << no_threads << " threads\n";
	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	perform_alignment_impl(clip_char, verbose);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	double duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
	std::cout << "   Alignment took "
			  << static_cast<uint32_t>(duration) / 3600 << "h "
			  << (static_cast<uint32_t>(duration) / 60) % 60 << "min "
			  << static_cast<uint32_t>(duration) % 60 << "s ("
			  << number_of_reads() / (duration * no_threads) << " reads/thread/s).\n";
}

template <typename T>
void single_end_aligner<T>::perform_alignment_impl(char clip_char, bool verbose) noexcept
{
	std::cout << "   Aligning reads from " << m_file_name << " (" << m_reads.size() << ")" << std::flush;
	alignment_impl(m_reads, m_profile_hmm, clip_char, verbose);
	std::cout << '\n';
}

template <typename T>
void paired_end_aligner<T>::perform_alignment_impl(char clip_char, bool verbose) noexcept
{
	single_end_aligner<T>::perform_alignment_impl(clip_char, verbose);

	std::cout << "   Aligning reads from " << m_file_name2 << " (" << m_reads2.size() << ")" << std::flush;
	alignment_impl(m_reads2, m_profile_hmm, clip_char, verbose);
	std::cout << '\n';
}

// write alignment
template <typename T>
void single_end_aligner<T>::write_alignment_to_file(const std::string& output_file_name, const std::string& rejects_file_name) noexcept
{
	std::cout << "6) Writing ";
	write_alignment_to_file_impl(output_file_name, rejects_file_name);
	std::cout << '\n';
}

void write_header(std::ostream& output, int profile_length, int argc, char** argv) noexcept
{
	output
		<< "@HD\tVN:1.5\tSO:queryname\n"
		<< "@SQ\tSN:CONSENSUS\tLN:" << profile_length << '\n'
		<< "@PG\tID:" PACKAGE_NAME "\tPN:" PACKAGE_NAME "\tVN:" PACKAGE_VERSION "\tCL:" << argv[0];

	for (int i = 1; i < argc; ++i)
	{
		output << ' ' << argv[i];
	}
	output << '\n';
}

template <typename T>
void single_end_aligner<T>::write_alignment_to_file_impl(const std::string& output_file_name, const std::string& rejects_file_name) noexcept
{
	std::cout << "single-end alignment\n";
	std::size_t proper_reads = 0;
	std::size_t too_short_reads = 0;

	std::ofstream output_file(output_file_name);
	std::ofstream rejects_file(rejects_file_name);

	write_header(output_file, m_profile_hmm.get_length_profile(), m_argc, m_argv);
	write_header(rejects_file, m_profile_hmm.get_length_profile(), m_argc, m_argv);

	for (const auto& i : m_reads)
	{
		((i.m_sam_record.mapped_length < m_min_mapped_length) ? (++too_short_reads, rejects_file) : (++proper_reads, output_file)) << i;
	}

	std::cout
		<< "\n\tToo short:               " << too_short_reads << '\n'
		<< '\n'
		<< "\tReads                    " << proper_reads << " (" << std::fixed << std::setprecision(1) << static_cast<double>(proper_reads) * 100 / m_reads.size() << "%)\n";
}

template <typename T>
void paired_end_aligner<T>::write_alignment_to_file_impl(const std::string& output_file_name, const std::string& rejects_file_name) noexcept
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
		if (*it1 < *it2)
		{
			++it1;
		}
		else
		{
			if (!(*it2 < *it1))
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
						if ((it2->m_sam_record.forward) && (it1->m_sam_record.forward == false))
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

				if (it1->m_sam_record.forward == false)
				{
					it2->m_sam_record.FLAG |= 0x20;
				}
				if (it2->m_sam_record.forward == false)
				{
					it1->m_sam_record.FLAG |= 0x20;
				}

				it1->m_sam_record.RNEXT = '=';
				it2->m_sam_record.RNEXT = '=';

				it1->m_sam_record.PNEXT = it2->m_sam_record.POS;
				it2->m_sam_record.PNEXT = it1->m_sam_record.POS;

				TLEN = static_cast<int32_t>(it2->m_sam_record.POS + it2->m_sam_record.segment_length) - static_cast<int32_t>(it1->m_sam_record.POS);
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

	write_header(output_file, m_profile_hmm.get_length_profile(), m_argc, m_argv);
	write_header(rejects_file, m_profile_hmm.get_length_profile(), m_argc, m_argv);

	it1 = m_reads.begin();
	it2 = m_reads2.begin();

	for (; (it1 != end1) && (it2 != end2);)
	{
		// check whether reads are too short
		if (it1->m_sam_record.mapped_length < m_min_mapped_length)
		{
			++too_short_reads;
			rejects_file << *it1;
			++it1;
			continue;
		}
		if (it2->m_sam_record.mapped_length < m_min_mapped_length)
		{
			++too_short_reads;
			rejects_file << *it2;
			++it2;
			continue;
		}

		// check whether reads are properly paired
		if (*it1 < *it2)
		{
			++unpaired_reads;
			unpaired_file << *it1;
			++it1;
		}
		else
		{
			if (!(*it2 < *it1))
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
		(it1->m_sam_record.mapped_length < m_min_mapped_length ? (++too_short_reads, rejects_file) : (++unpaired_reads, unpaired_file)) << *it1;
	}
	for (; it2 != end2; ++it2)
	{
		(it2->m_sam_record.mapped_length < m_min_mapped_length ? (++too_short_reads, rejects_file) : (++unpaired_reads, unpaired_file)) << *it2;
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

#endif /* ALIGNER_IMPL_HPP */
