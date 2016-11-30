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
public:
	// MARGINAL PROBABILITY of a sequence
	static T logLik(const reference_genome<T>& parameters, const boost::string_ref& sequence);

	static inline T Lik(const reference_genome<T>& parameters, const boost::string_ref& sequence)
	{
		return exp_base(logLik(parameters, sequence));
	}

	// OPTIMAL ALIGNMENT of a sequence
	static T viterbi(
		const reference_genome<T>& parameters,
		const boost::string_ref& sequence,
		uint32_t ref_start,
		uint32_t ref_end,
		std::vector<minimal_alignment>& alignments) noexcept;
};
} // unnamed namespace

#include "hmmalign_impl.hpp"

#endif /* NGSHMMALIGN_HMMALIGN_HPP */