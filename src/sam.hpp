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

/* sam_entry */
struct sam_entry
{
	// SAM fields
	uint16_t FLAG = 0;
	uint32_t POS = 0;
	uint16_t MAPQ = 254;
	std::string CIGAR;
	uint32_t PNEXT = 0;
	int32_t TLEN = 0;
	char RNEXT = '*';

	// custom fields
	bool forward = true;
	int64_t Score;
	std::size_t mapped_length;
	std::size_t segment_length;
	std::size_t Left_clip_length;
	std::size_t Right_clip_length;

	sam_entry() = default;
	sam_entry(const sam_entry&) = default;
	sam_entry(sam_entry&&) = default;
	sam_entry& operator=(const sam_entry& other) = default;
	sam_entry& operator=(sam_entry&& other) = default;

	sam_entry(uint32_t POS_, std::string&& CIGAR_, int64_t Score_, std::size_t mapped_length_, std::size_t segment_length_, std::size_t Left_clip_length_, std::size_t Right_clip_length_)
		: POS(POS_), CIGAR(std::move(CIGAR_)), Score(Score_), mapped_length(mapped_length_), segment_length(segment_length_), Left_clip_length(Left_clip_length_), Right_clip_length(Right_clip_length_)
	{
	}
};

#endif /* SAM_HPP */
