#ifndef NGSHMMALIGN_FASTQ_HPP
#define NGSHMMALIGN_FASTQ_HPP

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
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/iostreams/device/mapped_file.hpp>

#include "debug.hpp"

namespace
{

template <typename T>
void fastq_read(const std::string& input_file, std::vector<T>& reads, const bool paired_end, const bool second_in_pair) noexcept
{
	const boost::iostreams::mapped_file_source file(input_file);
	fastq_read(file, reads, paired_end, second_in_pair);
}

template <typename T>
void fastq_read(const boost::iostreams::mapped_file_source& input_file, std::vector<T>& reads, const bool paired_end, const bool second_in_pair) noexcept
{
	uint8_t phase = 0;
	std::size_t line = 0;

	const char* id_ptr = nullptr;
	std::size_t id_len = 0;

	const char* seq_ptr = nullptr;
	std::size_t seq_len = 0;

	const char* qual_ptr = nullptr;
	std::size_t qual_len = 0;

	bool found_space = false;
	bool read_second_in_pair = false;
	char first_letter_after_space = '\0';
	int header_length = 0;
	int spacer_length = 0;

	DEBUG_TRACE(std::cerr << std::boolalpha);

	const char* start_line_ptr = input_file.data();
	const char* const end_ptr = input_file.data() + input_file.size();
	for (const char* i = input_file.data(); i < end_ptr; ++i)
	{
		if (*i == '\n')
		{
			++line;

			if (i != start_line_ptr)
			{
				switch (phase)
				{
					case 0: // header line
						// @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
						//                                         ^
						//                identifies first(1) or second(2) in pair
						if (start_line_ptr[0] != '@')
						{
							std::cerr << "ERROR: Malformed FASTQ header line " << input_file << " on line " << line << " does not start with '@'." << std::endl;
							exit(EXIT_FAILURE);
						}

						id_ptr = start_line_ptr + 1;
						header_length = i - id_ptr;
						if (found_space == true)
						{
							// identifier contains space, can possibly
							// extract mate-pair information
							if ((paired_end == true) && (header_length >= id_len + 2))
							{
								// has at least one character overhanging
								first_letter_after_space = *(id_ptr + id_len + 1);
								read_second_in_pair = (first_letter_after_space == '2');
							}
							else
							{
								read_second_in_pair = (paired_end == true) && (second_in_pair == true);
							}
						}
						else
						{
							// identifier does not contain space,
							// hence does not allow for determining mate-pair information
							id_len = header_length;
							read_second_in_pair = (paired_end == true) && (second_in_pair == true);
						}

						DEBUG_TRACE(std::cerr << "Full Header:     " << boost::string_ref(start_line_ptr, header_length + 1) << " (len=" << header_length + 1 << ")" << std::endl);
						DEBUG_TRACE(std::cerr << "Header(I):       " << boost::string_ref(id_ptr, id_len) << " (len=" << id_len << ")" << std::endl);
						DEBUG_TRACE(std::cerr << "Contains space:  " << found_space << std::endl);
						DEBUG_TRACE(if (found_space) std::cerr
							<< "Header(II):      " << boost::string_ref(id_ptr + id_len + 1, header_length - id_len - 1) << " (len=" << header_length - id_len - 1 << ")" << std::endl
							<< "Letter:          " << (first_letter_after_space ? std::string(1, first_letter_after_space) : std::string("N/A")) << std::endl);
						DEBUG_TRACE(std::cerr << "Second in pair:  " << read_second_in_pair << std::endl);
						break;

					case 1: // DNA sequence
						seq_ptr = start_line_ptr;
						seq_len = i - seq_ptr;

						DEBUG_TRACE(std::cerr << "DNA sequence:    " << boost::string_ref(seq_ptr, seq_len) << " (len=" << seq_len << ")" << std::endl);
						break;

					case 2: // '+' with possibly id following again
						if (start_line_ptr[0] != '+')
						{
							std::cerr << "ERROR: Malformed FASTQ spacer line " << input_file << " on line " << line << " does not start with '+'." << std::endl;
							exit(EXIT_FAILURE);
						}

						spacer_length = i - start_line_ptr - 1;
						if ((spacer_length) && (spacer_length != header_length))
						{
							// should contain copy of sequence
							std::cerr << "ERROR: Malformed FASTQ spacer line " << input_file << " on line " << line << " does not have same length (" << spacer_length << ") as header length (" << header_length << ")." << std::endl;
							exit(EXIT_FAILURE);
						}

						DEBUG_TRACE(std::cerr << "Spacer line:     " << boost::string_ref(start_line_ptr + 1, spacer_length) << " (len=" << spacer_length << ")" << std::endl);
						break;

					case 3: // Phred scores
						qual_ptr = start_line_ptr;
						qual_len = i - qual_ptr;

						if (seq_len != qual_len)
						{
							std::cerr << "ERROR: Sequence and Phred qualities have different lengths on line " << line << '.' << std::endl;
							exit(EXIT_FAILURE);
						}

						reads.emplace_back(id_ptr, id_len, seq_ptr, seq_len, qual_ptr, qual_len, read_second_in_pair);
						DEBUG_TRACE(std::cerr << "Phred/Qual line: " << boost::string_ref(qual_ptr, qual_len) << " (len=" << qual_len << ")" << std::endl
											  << std::endl);

						// reset for next read
						found_space = false;
						first_letter_after_space = '\0';
						break;
				}

				++phase;
				phase %= 4;
			}

			start_line_ptr = i + 1;
		}
		else
		{
			const char cur_letter = *i;
			switch (phase)
			{
				case 0:
					// somewhere in the middle of the header line
					if ((cur_letter == ' ') && (found_space == false))
					{
						// found the first space in the header line,
						// need to clip off everything beyond
						found_space = true;
						id_len = i - start_line_ptr - 1;
					}
					break;

				case 1:
					switch (cur_letter)
					{
						case 'A':
						case 'C':
						case 'G':
						case 'T':
						case 'N':
							break;

						case 'a':
						case 'c':
						case 'g':
						case 't':
						case 'n':
							std::cerr << "ERROR: Malformed FASTQ sequence line " << input_file << " on line " << line << " contains lowercase DNA character '" << cur_letter << "'. DNA characters need to be uppercase." << std::endl;
							exit(EXIT_FAILURE);
							break;

						default:
							std::cerr << "ERROR: Malformed FASTQ sequence line " << input_file << " on line " << line << " contains invalid character '" << cur_letter << "'. Allowed DNA characters are 'A', 'C', 'G', 'T' and 'N'." << std::endl;
							exit(EXIT_FAILURE);
							break;
					}
					break;

				case 3:
					if ((cur_letter < 33) || (cur_letter > 126))
					{
						std::cerr << "ERROR: Malformed FASTQ quality line " << input_file << " on line " << line << " contains invalid quality character '" << cur_letter << "'. The range has to start at '!' (= 0) and end at '~' (= 93), inclusive." << std::endl;
						exit(EXIT_FAILURE);
					}

					break;
			}
		}
	}
}
} // unnamed namespace

#endif /* NGSHMMALIGN_FASTQ_HPP */
