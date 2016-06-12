#ifndef FASTQ_HPP
#define FASTQ_HPP

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
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/iostreams/device/mapped_file.hpp>

namespace
{

template <typename T>
std::vector<T> fastq_read(const std::string& input_file) noexcept
{
	std::vector<T> result;

	std::string id, seq, qual;

	std::string temp;
	std::ifstream input(input_file);
	uint64_t phase = 0, line = 0;

	auto inserter = [&result](std::string& id, std::string& seq, std::string& qual)
	{
		if (!id.empty())
		{
			result.emplace_back(id.substr(1, id.find_first_of(' ') - 1), std::move(seq), std::move(qual));
			id.clear();
			seq.clear();
			qual.clear();
		}
	};

	if (input.is_open())
	{
		while (input.good())
		{
			std::getline(input, temp);
			temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());
			++line;

			if (temp.empty())
				continue;

			switch (phase)
			{
				case 0: // header line
					inserter(id, seq, qual);
					id = std::move(temp);
					break;

				case 1: // DNA sequence
					seq = std::move(temp);
					break;

				case 2: // '+' with possibly id following again
					if (temp[0] != '+')
					{
						std::cerr << "ERROR: Malformed FASTQ file " << input_file << " on line " << line << ".\n";
						exit(EXIT_FAILURE);
					}
					break;

				case 3: // Phred scores
					if (temp.length() != seq.length())
					{
						std::cerr << "ERROR: Sequence and Phred qualities have different lengths on line " << line << ".\n";
						exit(EXIT_FAILURE);
					}
					qual = std::move(temp);
					break;
			}

			++phase;
			phase %= 4;
		}
	}
	else
	{
		std::cerr << "Input file '" << input_file << "' could not be opened!\n";
		exit(EXIT_FAILURE);
	}
	input.close();

	// add the last line
	if ((!id.empty()) && (!seq.empty()) && (!qual.empty()))
	{
		inserter(id, seq, qual);
	}

	return result;
}

template <typename T>
std::vector<T> fastq_read(const boost::iostreams::mapped_file_source& input_file) noexcept
{
	std::vector<T> result;

	uint8_t phase = 0;
	std::size_t line = 0;

	const char* id_ptr;
	std::size_t id_len;
	bool found_space = false;

	const char* seq_ptr;
	std::size_t seq_len;

	const char* qual_ptr;
	std::size_t qual_len;

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
						id_ptr = start_line_ptr + 1;
						if (found_space == false)
						{
							id_len = i - id_ptr;
						}
						if (start_line_ptr[0] != '@')
						{
							std::cerr << "ERROR: Malformed FASTQ header " << input_file << " on line " << line << ".\n";
							exit(EXIT_FAILURE);
						}
						break;

					case 1: // DNA sequence
						seq_ptr = start_line_ptr;
						seq_len = i - seq_ptr;
						break;

					case 2: // '+' with possibly id following again
						if (start_line_ptr[0] != '+')
						{
							std::cerr << "ERROR: Malformed FASTQ spacer line " << input_file << " on line " << line << ".\n";
							exit(EXIT_FAILURE);
						}
						break;

					case 3: // Phred scores
						qual_ptr = start_line_ptr;
						qual_len = i - qual_ptr;

						if (seq_len != qual_len)
						{
							std::cerr << "ERROR: Sequence and Phred qualities have different lengths on line " << line << ".\n";
							exit(EXIT_FAILURE);
						}

						found_space = false;
						result.emplace_back(id_ptr, id_len, seq_ptr, seq_len, qual_ptr, qual_len);
						break;
				}

				++phase;
				phase %= 4;
			}

			start_line_ptr = i + 1;
		}
		else
		{
			// somewhere in the middle of a line
			if ((*i == ' ') && (phase == 0) && (found_space == false))
			{
				// found the first space in the header line,
				// need to clip off everything beyond
				found_space = true;
				id_len = i - start_line_ptr - 1;
			}
		}
	}

	return result;
}
}

#endif /* FASTQ_HPP */