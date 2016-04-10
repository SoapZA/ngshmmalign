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
			result.emplace_back(id.substr(1, id.find_first_of(' ')), std::move(seq), std::move(qual));
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

#endif /* FASTQ_HPP */
