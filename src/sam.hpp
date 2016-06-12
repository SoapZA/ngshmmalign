#ifndef SAM_HPP
#define SAM_HPP

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
#include <algorithm>
#include <ostream>
#include <string>
#include <utility>

#include "utility_functions.hpp"

namespace
{

enum class clip_mode
{
	soft,
	hard,
	HARD
};

using CIGAR_vec = std::vector<std::pair<char, uint32_t>>;

/* sam_entry */
struct sam_entry
{
	// SAM fields
	uint16_t FLAG = 0;
	const char* RNAME = nullptr;
	uint32_t POS = 0;
	uint16_t MAPQ = 254;

	CIGAR_vec CIGAR;

	uint32_t PNEXT = 0;
	int32_t TLEN = 0;
	char RNEXT = '*';

	// auxiliary fields
	int64_t SCORE = lower_limit;
	uint32_t NM = 1000;

	// custom fields
	bool m_forward = true;

	uint32_t m_mapped_length = 0;
	uint32_t m_segment_length = 0;

	uint32_t m_left_clip_length = 0;
	uint32_t m_right_clip_length = 0;

	clip_mode m_clip = clip_mode::hard;

	sam_entry() = default;
	sam_entry(const sam_entry&) = default;
	sam_entry(sam_entry&&) = default;
	sam_entry& operator=(const sam_entry& other) = default;
	sam_entry& operator=(sam_entry&& other) = default;

	sam_entry(
		const char* RNAME_,
		uint32_t POS_,
		CIGAR_vec&& CIGAR_,
		int64_t SCORE_,
		uint32_t NM_,
		uint32_t mapped_length_,
		uint32_t segment_length_,
		uint32_t left_clip_length_,
		uint32_t right_clip_length_)
		: RNAME(RNAME_),
		  POS(POS_),
		  CIGAR(std::move(CIGAR_)),
		  SCORE(SCORE_),
		  NM(NM_),
		  m_mapped_length(mapped_length_),
		  m_segment_length(segment_length_),
		  m_left_clip_length(left_clip_length_),
		  m_right_clip_length(right_clip_length_)
	{
	}
};
}

#endif /* SAM_HPP */
