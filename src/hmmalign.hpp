#ifndef NGSHMMALIGN_HMMALIGN_HPP
#define NGSHMMALIGN_HMMALIGN_HPP

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

#include <cmath>
#include <list>
#include <random>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "dna_array.hpp"
#include "reference.hpp"

namespace
{

// first: operation
//
using CIGAR_vec = std::vector<std::pair<char, uint32_t>>;

struct minimal_alignment
{
	CIGAR_vec m_CIGAR;

	// 0-based position of first mapped base
	int32_t m_POS;

	// number of bases clipped on left
	int32_t m_left_clip_length = 0;

	// number of bases clipped on right
	int32_t m_right_clip_length = 0;

	// number of mapped bases
	// = (number of bases) - (number of clipped bases)
	// = (number of M ops) + (number of I ops)
	int32_t m_mapped_bases = 0;

	// length of segment on genome from first mapped base to last mapped base
	// = (number of M ops) + (number of D ops)
	int32_t m_segment_length = 0;

	minimal_alignment(const std::vector<char>& alignment_, uint32_t POS_) noexcept;
};

template <typename T>
class hmmalign
{
private:
	static_assert(std::is_integral<T>::value, "T needs to be an integral type!\n");

	// left struct
	struct left_clip_struct
	{
		T score;

		bool from_begin : 1;
		bool from_left_clip : 1;
	};

	// main matrix
	struct DP_entry_struct
	{
		struct match_struct
		{
			T score;

			bool from_begin : 1;
			bool from_left_clip : 1;
			bool from_match : 1;
			bool from_insertion : 1;
			bool from_deletion : 1;
		} match;

		struct insertion_struct
		{
			T score;

			bool from_match : 1;
			bool from_insertion : 1;
		} insertion;

		struct deletion_struct
		{
			T score;

			bool from_match : 1;
			bool from_deletion : 1;
		} deletion;
	};

	// right clips
	struct right_clip_struct
	{
		T score;

		bool from_right_clip : 1;
		boost::dynamic_bitset<> from_match;
	};

	// terminal
	struct end_struct
	{
		T score;

		bool from_right_clip : 1;
		boost::dynamic_bitset<> from_match;
	};

	using genome_index_type = typename std::vector<typename reference_genome<T>::template trans_matrix<T>>::size_type;
	using seq_index_type = std::string::size_type;

	// memory-related variables
	seq_index_type m_seq_max = 0;
	genome_index_type m_region_max = 0;

	std::unique_ptr<left_clip_struct[]> left_clip_matrix;
	std::unique_ptr<DP_entry_struct[]> DP_matrix;
	std::unique_ptr<right_clip_struct[]> right_clip_matrix;
	end_struct end;

	// current alignment data
	uint32_t reference_start;
	uint32_t reference_end;
	seq_index_type seq_L;
	genome_index_type region_L;

	inline void viterbi_memory_requirements(const int32_t seq_length, const int32_t region_length)
	{
		bool reallocate = false;

		if (seq_length > m_seq_max)
		{
			reallocate = true;
			m_seq_max = seq_length;
		}

		if (region_length > m_region_max)
		{
			reallocate = true;
			m_region_max = region_length;
		}

		if (reallocate)
		{
			left_clip_matrix.reset(new left_clip_struct[m_seq_max - 1]);
			DP_matrix.reset(new DP_entry_struct[m_seq_max * m_region_max]);
			right_clip_matrix.reset(new right_clip_struct[m_seq_max - 1]);
		}
	}

public:
	// MARGINAL PROBABILITY of a sequence
	static T logLik(const reference_genome<T>& parameters, const boost::string_ref& sequence);

	static inline T Lik(const reference_genome<T>& parameters, const boost::string_ref& sequence)
	{
		return exp_base(logLik(parameters, sequence));
	}

	hmmalign() = default;

	// OPTIMAL ALIGNMENT of a sequence
	T viterbi(
		const reference_genome<T>& parameters,
		const boost::string_ref& sequence,
		const uint32_t ref_start,
		const uint32_t ref_end,
		std::vector<minimal_alignment>& alignments) noexcept;

private:
	inline void viterbi_initialize_impl(
		const reference_genome<T>& parameters,
		const boost::string_ref& sequence) noexcept;
	inline void viterbi_recursion_impl(
		const reference_genome<T>& parameters,
		const boost::string_ref& sequence) noexcept;
	inline T viterbi_backtrack_impl(
		const reference_genome<T>& parameters,
		const boost::string_ref& sequence,
		std::vector<minimal_alignment>& alignments) noexcept;
};
} // unnamed namespace

#include "hmmalign_impl.hpp"

#endif /* NGSHMMALIGN_HMMALIGN_HPP */
