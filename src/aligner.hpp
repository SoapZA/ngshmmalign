#ifndef ALIGNER_HPP
#define ALIGNER_HPP

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

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

template <typename T, std::size_t N>
class dna_array;
template <typename T>
struct background_rates;

#include "hmmalign.hpp"
#include "sam.hpp"

extern int no_threads;

struct fastq_entry
{
	std::string m_id;
	std::string m_seq;
	std::string m_qual;

	// ctor:
	template <typename T1, typename T2, typename T3>
	fastq_entry(T1&& id_, T2&& seq_, T3&& qual_);

	fastq_entry() = default;
	fastq_entry(const fastq_entry& other) = default;
	fastq_entry(fastq_entry&& other) = default;
	fastq_entry& operator=(const fastq_entry& other) = default;
	fastq_entry& operator=(fastq_entry&& other) = default;
};

struct read_entry
{
	// ctor:
	template <typename T1, typename T2, typename T3>
	read_entry(T1&& id_, T2&& seq_, T3&& qual_);

	read_entry() = default;
	read_entry(const read_entry& other) = default;
	read_entry(read_entry&& other) = default;
	read_entry& operator=(const read_entry& other) = default;
	read_entry& operator=(read_entry&& other) = default;

	inline bool operator<(const read_entry& other_read_entry) const noexcept;

	// output
	friend std::ostream& operator<<(std::ostream&, const read_entry&) noexcept;

	// members:
	fastq_entry m_fastq_record;
	sam_entry m_sam_record;
};

template <typename T>
class single_end_aligner
{
public:
	// ctor
	single_end_aligner(const std::string& file_name_, uint32_t min_mapped_length_, int argc_, char** argv_) noexcept;

	single_end_aligner() = delete;
	single_end_aligner(const single_end_aligner& other) = delete;
	single_end_aligner(single_end_aligner&& other) = delete;
	single_end_aligner& operator=(const single_end_aligner& other) = delete;
	single_end_aligner& operator=(single_end_aligner&& other) = delete;

	// dtor
	virtual ~single_end_aligner() = default;

	// mutators
	template <typename V>
	void set_parameters(
		const std::vector<dna_array<V, 5>>& allel_freq_,
		const std::vector<V>& vec_M_to_D_p_,
		const std::vector<V>& vec_D_to_D_p_,
		background_rates<V> error_rates) noexcept;

	// reads
	void sort_reads() noexcept;

	// alignment
	void perform_alignment(char clip_char, bool verbose) noexcept;

	// write alignment to output
	// output
	void write_alignment_to_file(const std::string& output_file_name, const std::string& rejects_file_name) noexcept;

protected:
	int m_argc;
	char** m_argv;

	uint32_t m_min_mapped_length;
	std::string m_file_name;
	std::vector<read_entry> m_reads;
	hmmalign<T> m_profile_hmm;

	// helper functions
	std::size_t number_of_reads() const noexcept;
	uint32_t get_length_profile() const noexcept;

	// reads
	virtual void sort_reads_impl() noexcept;
	virtual std::size_t number_of_reads_impl() const noexcept;

	// alignment
	virtual void perform_alignment_impl(char clip_char, bool verbose) noexcept;

	// write alignment to output
	virtual void write_alignment_to_file_impl(const std::string& output_file_name, const std::string& rejects_file_name) noexcept;
};

template <typename T>
class paired_end_aligner : public single_end_aligner<T>
{
public:
	// ctor
	paired_end_aligner(const std::string& file_name_, const std::string& file_name2_, bool write_unpaired, uint32_t min_mapped_length_, int argc_, char** argv_);

	paired_end_aligner() = delete;
	paired_end_aligner(const paired_end_aligner& other) = delete;
	paired_end_aligner(paired_end_aligner&& other) = delete;
	paired_end_aligner& operator=(const paired_end_aligner& other) = delete;
	paired_end_aligner& operator=(paired_end_aligner&& other) = delete;

private:
	using single_end_aligner<T>::m_argc;
	using single_end_aligner<T>::m_argv;
	using single_end_aligner<T>::m_min_mapped_length;
	using single_end_aligner<T>::m_file_name;
	using single_end_aligner<T>::m_reads;
	using single_end_aligner<T>::m_profile_hmm;

	std::string m_file_name2;
	std::vector<read_entry> m_reads2;
	bool m_write_unpaired;

	// reads
	virtual void sort_reads_impl() noexcept override;
	virtual std::size_t number_of_reads_impl() const noexcept override;

	// alignment
	virtual void perform_alignment_impl(char clip_char, bool verbose) noexcept override;

	// write alignment to output
	virtual void write_alignment_to_file_impl(const std::string& output_file_name, const std::string& rejects_file_name) noexcept override;
};

#include "aligner_impl.hpp"

#endif /* ALIGNER_HPP */
