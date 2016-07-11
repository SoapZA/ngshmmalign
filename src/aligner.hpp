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

#include <boost/utility/string_ref.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include "reference.hpp"
#include "hmmalign.hpp"
#include "sam.hpp"

extern int num_threads;

namespace
{

struct fastq_entry
{
	using string_rep = boost::string_ref;

	string_rep m_id;
	string_rep m_seq;
	string_rep m_qual;

	fastq_entry() = default;
	fastq_entry(const fastq_entry& other) = default;
	fastq_entry(fastq_entry&& other) = default;
	fastq_entry& operator=(const fastq_entry& other) = default;
	fastq_entry& operator=(fastq_entry&& other) = default;

	fastq_entry(
		const char* id_ptr,
		std::size_t id_len,
		const char* seq_ptr,
		std::size_t seq_len,
		const char* qual_ptr,
		std::size_t qual_len);
};

struct read_entry
{
	// ctor:
	read_entry(
		const char* id_ptr,
		std::size_t id_len,
		const char* seq_ptr,
		std::size_t seq_len,
		const char* qual_ptr,
		std::size_t qual_len);

	read_entry() = default;
	read_entry(const read_entry& other) = default;
	read_entry(read_entry&& other) = default;
	read_entry& operator=(const read_entry& other) = default;
	read_entry& operator=(read_entry&& other) = default;

	// output
	friend std::ostream& operator<<(std::ostream&, const read_entry&) noexcept;

	void print_sequence(std::ostream& output, bool clip_bases) const noexcept;
	void print_qual(std::ostream& output, bool clip_bases) const noexcept;

	// members:
	fastq_entry m_fastq_record;
	sam_entry m_sam_record;
};

struct partioned_genome;

template <typename T>
class single_end_aligner
{
public:
	static std::unique_ptr<single_end_aligner<T>> create_aligner_instance(const std::vector<std::string>& input_files, int32_t min_mapped_length, int argc, const char** argv, bool write_unpaired) noexcept;

	// 1. ctor
	single_end_aligner(int32_t min_mapped_length_, int argc_, const char** argv_) noexcept;

	single_end_aligner() = delete;
	single_end_aligner(const single_end_aligner& other) = delete;
	single_end_aligner(single_end_aligner&& other) = delete;
	single_end_aligner& operator=(const single_end_aligner& other) = delete;
	single_end_aligner& operator=(single_end_aligner&& other) = delete;

	// 2. load reads
	void load_reads(const std::vector<std::string>& input_files) noexcept;

	// 3. load parameters
	void load_parameters(const std::string& input_file, background_rates& error_rates, bool ambig_bases_unequal_weight) noexcept;

	// 4. sort reads
	void sort_reads() noexcept;

	// 5. perform parameter estimation
	void estimate_parameters(const std::string& data_root, const std::string& mafft, const background_rates& error_rates, uint64_t seed, bool verbose, bool keep_mafft_files, bool ambig_bases_unequal_weight) noexcept;

	// 6. perform alignment
	void perform_alignment(const std::string& reference_genome_name, clip_mode clip, uint64_t seed, bool exhaustive, bool verbose, bool differentiate_match_state) noexcept;

	// 7. write alignment to output
	void write_alignment_to_file(const std::string& data_root, const std::string& output_file_name, const std::string& rejects_file_name) noexcept;

	// 8. dtor
	virtual ~single_end_aligner() = default;

protected:
	// 2. load reads
	virtual void load_reads_impl(const std::vector<std::string>& input_files) noexcept;

	// 3. load parameters
	uint32_t get_length_profile() const noexcept;

	// 4. sort reads
	virtual void sort_reads_impl() noexcept;

	// 5. perform parameter estimation
	virtual void estimate_parameters_impl(partioned_genome& separated_reads, std::default_random_engine& generator) noexcept;

	// 6. perform alignment
	virtual std::size_t number_of_reads() const noexcept;
	virtual void perform_alignment_impl(clip_mode clip, uint64_t seed, bool exhaustive, bool verbose, bool differentiate_match_state) noexcept;

	// 7. write alignment to output
	virtual void write_alignment_to_file_impl(const std::string& data_root, const std::string& output_file_name, const std::string& rejects_file_name) noexcept;

	// command line parameters
	int m_argc;
	const char** m_argv;
	int m_phase = 0;

	// read data
	int32_t m_min_mapped_length;
	std::string m_read_file_name;
	boost::iostreams::mapped_file_source m_reads_mmap;
	std::vector<read_entry> m_reads;

	// profile HMM relevant members
	reference_genome<T> m_parameters;
	bool serialize = false;
};

template <typename T>
class paired_end_aligner : public single_end_aligner<T>
{
public:
	// 1. ctor
	paired_end_aligner(int32_t min_mapped_length_, int argc_, const char** argv_, bool write_unpaired) noexcept;

	paired_end_aligner() = delete;
	paired_end_aligner(const paired_end_aligner& other) = delete;
	paired_end_aligner(paired_end_aligner&& other) = delete;
	paired_end_aligner& operator=(const paired_end_aligner& other) = delete;
	paired_end_aligner& operator=(paired_end_aligner&& other) = delete;

private:
	// 2. load reads
	virtual void load_reads_impl(const std::vector<std::string>& input_files) noexcept override;

	// 4. sort reads
	virtual void sort_reads_impl() noexcept override;

	// 5. perform parameter estimation
	virtual void estimate_parameters_impl(partioned_genome& separated_reads, std::default_random_engine& generator) noexcept override;

	// 6. perform alignment
	virtual std::size_t number_of_reads() const noexcept override;
	virtual void perform_alignment_impl(clip_mode clip, uint64_t seed, bool exhaustive, bool verbose, bool differentiate_match_state) noexcept override;

	// 7. write alignment to output
	virtual void write_alignment_to_file_impl(const std::string& data_root, const std::string& output_file_name, const std::string& rejects_file_name) noexcept override;

public:
	// 8. dtor
	virtual ~paired_end_aligner() = default;

private:
	// read data
	std::string m_read_file_name2;
	boost::iostreams::mapped_file_source m_reads2_mmap;
	std::vector<read_entry> m_reads2;
	bool m_write_unpaired;

	using single_end_aligner<T>::m_argc;
	using single_end_aligner<T>::m_argv;

	using single_end_aligner<T>::m_min_mapped_length;
	using single_end_aligner<T>::m_read_file_name;
	using single_end_aligner<T>::m_reads;

	using single_end_aligner<T>::m_parameters;
};
}

#include "aligner_impl.hpp"

#endif /* ALIGNER_HPP */