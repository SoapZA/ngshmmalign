#ifndef REFERENCE_HPP
#define REFERENCE_HPP

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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "dna_array.hpp"

const std::map<const std::string, const char> ambig_to_wobble_base = {
	{ "", 'N' },
	{ "A", 'A' },
	{ "C", 'C' },
	{ "G", 'G' },
	{ "T", 'T' },
	{ "AC", 'M' },
	{ "AG", 'R' },
	{ "AT", 'W' },
	{ "CG", 'S' },
	{ "CT", 'Y' },
	{ "GT", 'K' },
	{ "ACG", 'V' },
	{ "ACT", 'H' },
	{ "AGT", 'D' },
	{ "CGT", 'B' },
	{ "ACGT", 'N' }
};

struct reference_haplotype
{
	std::string name;
	std::string sequence;

	// we use half-open intervals
	// i.e., for the enclosed-interval [start, end)
	std::string::size_type start;
	std::string::size_type end;
	double count;

	reference_haplotype(const std::string& id, std::string&& seq) noexcept : sequence(std::move(seq))
	{
		std::vector<std::string> split_vec;
		boost::split(split_vec, id, boost::is_any_of("_"), boost::token_compress_on);

		if (split_vec.size() == 1)
		{
			count = 1;
		}
		else
		{
			try
			{
				count = boost::lexical_cast<double>(split_vec[1]);
			}
			catch (boost::bad_lexical_cast&)
			{
				std::cerr << "ERROR: Count argument '" << split_vec[1] << "' is not an integral/floating point value! Aborting.\n";
				exit(EXIT_FAILURE);
			}
		}

		name = std::move(split_vec[0]);

		// find start
		start = 0;
		while (sequence[start] == '-')
		{
			++start;
		};

		// find end
		end = sequence.length() - 1;
		while (sequence[end] == '-')
		{
			--end;
		};
		++end;
	}
};

std::tuple<std::vector<dna_array<double, 5>>, std::vector<double>, std::vector<double>> build_parameters_from_msa(const std::vector<reference_haplotype>& refs, double ambig_threshold)
{
	const std::string::size_type L = refs[0].sequence.length();
	const std::vector<reference_haplotype>::size_type num_haps = refs.size();

	std::tuple<std::vector<dna_array<double, 5>>, std::vector<double>, std::vector<double>> result = std::make_tuple(
		std::vector<dna_array<double, 5>>(L, { 0.0, 0.0, 0.0, 0.0, 0.0 }),
		std::vector<double>(L - 1, 0),
		std::vector<double>(L - 1, 0));
	std::vector<dna_array<double, 5>>& E_p = std::get<0>(result);
	std::vector<double>& M_D_p = std::get<1>(result);
	std::vector<double>& D_D_p = std::get<2>(result);

	// 1.) check that all haplotypes have the same length
	for (std::vector<reference_haplotype>::size_type i = 1; i < num_haps; ++i)
	{
		if (refs[i].sequence.length() != L)
		{
			std::cerr << "ERROR: Haplotype '" << refs[i].name << "' does not have length L = " << L << "! Aborting.\n";
			exit(EXIT_FAILURE);
		}

		if (refs[i].count <= 0)
		{
			std::cerr << "ERROR: Haplotype '" << refs[i].name << "' has non-positive count " << refs[i].count << "! Aborting.\n";
			exit(EXIT_FAILURE);
		}
	}

	// 2.) check that all haplotypes have proper DNA letters
	char cur_base;
	for (std::vector<reference_haplotype>::size_type i = 0; i < num_haps; ++i)
	{
		for (std::string::size_type j = 0; j < L; ++j)
		{
			cur_base = refs[i].sequence[j];
			if ((cur_base != 'A') && (cur_base != 'C') && (cur_base != 'G') && (cur_base != 'T') && (cur_base != '-'))
			{
				std::cerr << "ERROR: Unknown base '" << cur_base << "' at position '" << j << "' of '" << refs[i].name << "'.\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	// 3.) start producing emission and transition tables
	double sum;
	double MD, sumM;
	double DD, sumD;
	bool onlyGap;

	std::string majority_bases;
	std::string majority_ref;
	double max_freq;
	std::default_random_engine rng;
	std::uniform_int_distribution<int> dist;
	char majority_base;

	std::string ambig_ref;
	std::string ambig_bases;

	for (std::string::size_type j = 0; j < L; ++j)
	{
		sum = 0;

		MD = 0;
		sumM = 0;
		DD = 0;
		sumD = 0;

		onlyGap = true;

		// first loop
		for (std::vector<reference_haplotype>::size_type i = 0; i < num_haps; ++i)
		{
			if ((refs[i].start <= j) && (j < refs[i].end))
			{
				cur_base = refs[i].sequence[j];
				onlyGap &= (cur_base == '-');

				// allele frequencies
				if (cur_base != '-')
				{
					sum += refs[i].count;
					E_p[j][cur_base] += refs[i].count;
				}

				// transition tables
				if (j < L - 1)
				{
					if (cur_base == '-')
					{
						// D->
						sumD += refs[i].count;
						if (refs[i].sequence[j + 1] == '-')
						{
							// ->D
							DD += refs[i].count;
						}
					}
					else
					{
						// M->
						sumM += refs[i].count;
						if (refs[i].sequence[j + 1] == '-')
						{
							// ->D
							MD += refs[i].count;
						}
					}
				}
			}
		}

		if (onlyGap)
		{
			std::cerr << "ERROR: Position '" << j << "' only has gaps.\n";
			exit(EXIT_FAILURE);
		}

		majority_bases.clear();
		max_freq = -1;
		ambig_bases.clear();
		// renormalize counts
		for (char k : { 'A', 'C', 'G', 'T' })
		{
			E_p[j][k] /= sum;

			// find the majority base
			if (E_p[j][k] >= max_freq)
			{
				if (E_p[j][k] > max_freq)
				{
					majority_bases.clear();
					max_freq = E_p[j][k];
				}

				majority_bases.push_back(k);
			}

			// find ambiguous base
			if (E_p[j][k] >= ambig_threshold)
			{
				ambig_bases.push_back(k);
			}
		}

		dist.param(std::uniform_int_distribution<int>::param_type(0, majority_bases.length() - 1));
		majority_base = majority_bases[dist(rng)];
		majority_ref.push_back(majority_base);

		ambig_ref.push_back(ambig_to_wobble_base.find(ambig_bases)->second);

		if (j < L - 1)
		{
			M_D_p[j] = MD / (sumM ? sumM : 1);
			D_D_p[j] = DD / (sumD ? sumD : 1);
		}
	}

	if (num_haps > 1)
	{
		std::ofstream output;
		output.open("ref_majority.fasta");
		output << ">CONSENSUS" << '\n' << majority_ref << '\n';
		output.close();

		output.open("ref_ambig.fasta");
		output << ">CONSENSUS" << '\n' << ambig_ref << '\n';
		output.close();
	}

	return result;
}

#endif /* REFERENCE_HPP */