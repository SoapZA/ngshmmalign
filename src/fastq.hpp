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
#include <boost/regex.hpp>

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
	std::size_t line = 1;

	const char* id_ptr = nullptr;
	std::size_t id_len = 0;

	const char* seq_ptr = nullptr;
	std::size_t seq_len = 0;

	const char* qual_ptr = nullptr;
	std::size_t qual_len = 0;

	bool read_second_in_pair = false;
	int header_length = 0;
	int spacer_length = 0;

	DEBUG_TRACE(std::cerr << std::boolalpha);

	// Identify Illumina headers
	// https://en.wikipedia.org/wiki/FASTQ_format
	//
	// 1.) Pre-1.8 format:
	//     @HWUSI-EAS100R:6:73:941:1973#0/1
	//
	//       HWUSI-EAS100R:6:73:941:1973#0 -> readname
	//       /1                            -> matepair
	//
	// 2.) 1.8+ format:
	//     @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	//
	//       EAS139:136:FC706VJ:2:2104:15343:197393 -> readname
	//       1                                      -> matepair
	//       :Y:18:ATCACG                           -> remaining auxiliary info (discarded)

	boost::regex illumina_old{ "([[:graph:]]*)/([[:digit:]])" };
	boost::regex illumina_new{ "([[:graph:]]*) ([[:digit:]]):[[:graph:]]*:[[:graph:]]*:[[:graph:]]*" };
	boost::cmatch header_matches;
	bool found_match = false;

	enum class illumina_format
	{
		OLD,
		NEW,
		UNKNOWN
	} header_format;

	const char* start_line_ptr = input_file.data();
	const char* const end_ptr = input_file.data() + input_file.size();
	for (const char* i = input_file.data(); i < end_ptr; ++i)
	{
		if (*i == '\n')
		{
			if (i != start_line_ptr)
			{
				switch (phase)
				{
					case 0: // header line
						if (*start_line_ptr != '@')
						{
							std::cerr << "ERROR: Malformed FASTQ header: line " << line << " does not start with '@'." << std::endl;
							exit(EXIT_FAILURE);
						}

						id_ptr = start_line_ptr + 1;
						header_length = i - id_ptr;

						found_match = ((header_format = illumina_format::NEW, boost::regex_search(id_ptr, i, header_matches, illumina_new)) || (header_format = illumina_format::OLD, boost::regex_search(id_ptr, i, header_matches, illumina_old)) || (header_format = illumina_format::UNKNOWN, false));

						read_second_in_pair = paired_end;
						if (found_match)
						{
							id_ptr = header_matches[1].first;
							id_len = header_matches[1].second - id_ptr;

							const char mate_pair_identifier = *(header_matches[2].first);
							switch (mate_pair_identifier)
							{
								case '1':
									read_second_in_pair = false;
									break;
								case '2':
									read_second_in_pair &= true;
									break;
								default:
									std::cerr << "ERROR: Malformed FASTQ header: line " << line << " contained invalid matepair specifier '" << mate_pair_identifier << "'." << std::endl;
									exit(EXIT_FAILURE);
									break;
							}
						}
						else
						{
							id_len = i - id_ptr;
							read_second_in_pair &= second_in_pair;
						}

						DEBUG_TRACE(std::cerr << "Full Header:     '" << boost::string_ref(id_ptr, header_length) << "' (len=" << header_length << ")" << std::endl);
						DEBUG_TRACE(std::cerr << "Format:           " << (header_format == illumina_format::NEW ? "New Illumina 1.8+" : (header_format == illumina_format::OLD ? "Old Illumina <1.8" : "Unknown format")) << std::endl);
						DEBUG_TRACE(std::cerr << "Readname:        '" << boost::string_ref(id_ptr, id_len) << "' (len=" << id_len << ")" << std::endl);
						DEBUG_TRACE(std::cerr << "Second in pair:   " << read_second_in_pair << std::endl);
						break;

					case 1: // DNA sequence
						seq_ptr = start_line_ptr;
						seq_len = i - seq_ptr;

						DEBUG_TRACE(std::cerr << "DNA sequence:    '" << boost::string_ref(seq_ptr, seq_len) << "' (len=" << seq_len << ")" << std::endl);
						break;

					case 2: // '+' with possibly id following again
						if (start_line_ptr[0] != '+')
						{
							std::cerr << "ERROR: Malformed FASTQ spacer: line " << line << " does not start with '+'." << std::endl;
							exit(EXIT_FAILURE);
						}

						spacer_length = i - start_line_ptr - 1;
						if ((spacer_length) && (spacer_length != header_length))
						{
							// should contain copy of sequence
							std::cerr << "ERROR: Malformed FASTQ spacer: line " << line << " does not have same length (" << spacer_length << ") as header length (" << header_length << ")." << std::endl;
							exit(EXIT_FAILURE);
						}

						DEBUG_TRACE(std::cerr << "Spacer line:     '" << boost::string_ref(start_line_ptr + 1, spacer_length) << "' (len=" << spacer_length << ")" << std::endl);
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
						DEBUG_TRACE(std::cerr << "Phred/Qual line: '" << boost::string_ref(qual_ptr, qual_len) << "' (len=" << qual_len << ")" << std::endl
											  << std::endl);
						break;
				}

				++phase;
				phase %= 4;
			}

			start_line_ptr = i + 1;
			++line;
		}
		else
		{
			const char cur_letter = *i;
			switch (phase)
			{
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
							std::cerr << "ERROR: Malformed FASTQ sequence: line " << line << " contains lowercase DNA character '" << cur_letter << "'. DNA characters need to be uppercase." << std::endl;
							exit(EXIT_FAILURE);
							break;

						default:
							std::cerr << "ERROR: Malformed FASTQ sequence: line " << line << " contains invalid character '" << cur_letter << "'. Allowed DNA characters are 'A', 'C', 'G', 'T' and 'N'." << std::endl;
							exit(EXIT_FAILURE);
							break;
					}
					break;

				case 3:
					if ((cur_letter < 33) || (cur_letter > 126))
					{
						std::cerr << "ERROR: Malformed FASTQ qualities: line " << line << " contains invalid quality character '" << cur_letter << "'. The range has to start at '!' (= 0) and end at '~' (= 93), inclusive." << std::endl;
						exit(EXIT_FAILURE);
					}

					break;
			}
		}
	}
}
} // unnamed namespace

#endif /* NGSHMMALIGN_FASTQ_HPP */
