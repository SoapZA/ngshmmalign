#ifndef NGSHMMALIGN_ALIGNER_HPP
#define NGSHMMALIGN_ALIGNER_HPP

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

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/utility/string_ref.hpp>

#include "hmmalign.hpp"
#include "reference.hpp"

extern int num_threads;

namespace
{

enum class clip_mode
{
	soft,
	hard,
	HARD
};

struct fastq_entry
{
	std::string m_id;
	std::string m_seq;
	std::string m_qual;
	bool m_second_in_pair;

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
		std::size_t qual_len,
		bool second_in_pair);
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
		std::size_t qual_len,
		bool second_in_pair);

	read_entry() = default;
	read_entry(const read_entry& other) = default;
	read_entry(read_entry&& other) = default;
	read_entry& operator=(const read_entry& other) = default;
	read_entry& operator=(read_entry&& other) = default;

	// output
	friend std::ostream& operator<<(std::ostream&, const read_entry&) noexcept;
	inline void print_sequence_qual(std::ostream& output, bool clip_bases, boost::string_ref str) const noexcept;

	// set and get optimal alignment
	inline void set_opt_aln(std::size_t i) noexcept;
	inline minimal_alignment& get_opt_aln() noexcept;
	inline const minimal_alignment& get_opt_aln() const noexcept;

	static clip_mode m_clip;
	// members:
	fastq_entry m_fastq_record;
	std::vector<minimal_alignment> m_cooptimal_alignments;

	// read status
	bool m_invalid = false;

	// PE - 0x1
	bool m_paired_read = false;
	// PE - 0x2
	bool m_mapped_in_proper_pair = false;
	// PE - 0x8
	//      0x10
	bool m_reverse_compl = false;
	// PE - 0x20
	bool m_mate_reverse_compl = false;

	// SAM fields (for optimal alignment)
	static boost::string_ref RNAME;
	static constexpr uint16_t MAPQ = 254;

	int32_t PNEXT = 0;
	int32_t TLEN = 0;
	char RNEXT = '*';

	// auxiliary fields
	uint32_t NM = 0;
	int64_t SCORE = lower_limit;
	std::string MD_tag;
};

struct partioned_genome;

template <typename T>
class single_end_aligner
{
public:
	static std::unique_ptr<single_end_aligner<T>> create_aligner_instance(const std::vector<std::string>& input_files, int32_t min_required_mapped_bases, int argc, const char** argv) noexcept;

	// 1. ctor + dtor
	single_end_aligner(int32_t min_required_mapped_bases_, int argc_, const char** argv_) noexcept;

	single_end_aligner() = delete;
	single_end_aligner(const single_end_aligner& other) = delete;
	single_end_aligner(single_end_aligner&& other) = delete;
	single_end_aligner& operator=(const single_end_aligner& other) = delete;
	single_end_aligner& operator=(single_end_aligner&& other) = delete;

	virtual ~single_end_aligner() = default;

	// 2. load reads
	void load_reads(const std::vector<std::string>& input_files) noexcept;

	// 3. load parameters
	void load_parameters(const std::string& input_file, background_rates& error_rates, bool ambig_bases_unequal_weight) noexcept;

	// 4. perform parameter estimation
	void estimate_parameters(const std::string& data_root, const std::string& mafft, const background_rates& error_rates, uint64_t seed, bool verbose, bool keep_mafft_files, bool ambig_bases_unequal_weight, int32_t num_threads, int32_t chunk_size) noexcept;

	// 5. perform alignment
	void perform_alignment(bool exhaustive, bool verbose, int32_t num_threads, int32_t chunk_size) noexcept;

	// 6. post-alignment processing
	void post_alignment_processing(bool differentiate_match_state, uint64_t seed, double min_freq, double error_rate, bool ambig_bases_unequal_weight) noexcept;

	// 7. write alignment to output
	void write_alignment_to_file(const clip_mode clip, const std::string& consensus_name, const std::string& data_root, const std::string& output_file_name, const std::string& rejects_file_name) noexcept;

protected:
	// 2. load reads
	virtual void load_reads_impl(const std::vector<std::string>& input_files) noexcept;

	// 3. load parameters
	uint32_t get_length_profile() const noexcept;

	// 5. perform alignment
	void perform_alignment_impl(bool exhaustive, bool verbose, int32_t num_threads, int32_t chunk_size) noexcept;

	// 6. post-alignment processing
	virtual void flag_reads() noexcept;

	// command line parameters
	int m_argc;
	const char** m_argv;
	int m_phase = 0;

	// read data
	int32_t m_min_required_mapped_bases;
	std::string m_read_file_name;

	std::vector<read_entry> m_reads;
	std::vector<read_entry*> m_good_reads;
	std::vector<read_entry*> m_bad_reads;

	// profile HMM relevant members
	reference_genome<T> m_parameters;
};

template <typename T>
class paired_end_aligner : public single_end_aligner<T>
{
public:
	// 1. ctor + dtor
	paired_end_aligner(int32_t min_required_mapped_bases_, int argc_, const char** argv_) noexcept;

	paired_end_aligner() = delete;
	paired_end_aligner(const paired_end_aligner& other) = delete;
	paired_end_aligner(paired_end_aligner&& other) = delete;
	paired_end_aligner& operator=(const paired_end_aligner& other) = delete;
	paired_end_aligner& operator=(paired_end_aligner&& other) = delete;

	virtual ~paired_end_aligner() = default;

protected:
	// 2. load reads
	virtual void load_reads_impl(const std::vector<std::string>& input_files) noexcept override;

	// 6. post-alignment processing
	virtual void flag_reads() noexcept override;

	// read data
	std::string m_read_file_name2;

	using single_end_aligner<T>::m_argc;
	using single_end_aligner<T>::m_argv;
	using single_end_aligner<T>::m_phase;

	using single_end_aligner<T>::m_min_required_mapped_bases;
	using single_end_aligner<T>::m_read_file_name;

	using single_end_aligner<T>::m_reads;
	using single_end_aligner<T>::m_good_reads;
	using single_end_aligner<T>::m_bad_reads;

	using single_end_aligner<T>::m_parameters;
};
} // unnamed namespace

#include "aligner_impl.hpp"

#endif /* NGSHMMALIGN_ALIGNER_HPP */
