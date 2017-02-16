#ifndef NGSHMMALIGN_HMMALIGN_IMPL_HPP
#define NGSHMMALIGN_HMMALIGN_IMPL_HPP

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

namespace
{

minimal_alignment::minimal_alignment(const std::vector<char>& alignment_, uint32_t POS_) noexcept
	: m_POS(POS_)
{
	// fully expanded alignment operation string
	// e.g.
	//       0   4
	// Ref:  AGCATCA---GCTAGCTAGCTGACT
	// Read:     TCACCCGCT---TAG
	//
	// yields the full-length CIGAR
	//       MMMIIIMMMDDDMMM
	// which eventually gets turned into
	//       3M3I3M3D3M
	// here segment_length = 12, i.e. end_POS = POS + segment_length
	//

	m_CIGAR.reserve(5);

	char op = alignment_.back();
	uint32_t count = 0;

	const auto it_end = alignment_.crend();
	for (auto it = alignment_.crbegin(); it != it_end; ++it)
	{
		if (*it != op)
		{
			switch (op)
			{
				case 'L':
					m_left_clip_length = count;
					break;

				case 'R':
					m_right_clip_length = count;
					break;

				case 'M':
					m_mapped_bases += count;

				case 'D':
					m_segment_length += count;
					m_CIGAR.emplace_back(op, count);
					break;

				case 'I':
					m_mapped_bases += count;
					m_CIGAR.emplace_back(op, count);
					break;

				default:
					break;
			}
			op = *it;
			count = 1;
		}
		else
		{
			++count;
		}
	}
}

template <typename T>
T hmmalign<T>::logLik(
	const reference_genome<T>& parameters,
	const boost::string_ref& sequence)
{
	static_assert(std::is_floating_point<T>::value, "T needs to be a floating point type!\n");
	using genome_index_type = typename std::vector<typename reference_genome<T>::template trans_matrix<T>>::size_type;
	using seq_index_type = std::string::size_type;

	const auto seq_L = sequence.length();

	// DP matrices
	struct DP_entry
	{
		T ln_match;
		T ln_insertion;
		T ln_deletion;
	};
	std::unique_ptr<T[]> left_clip_matrix(new T[seq_L - 1]);
	std::unique_ptr<DP_entry[]> DP_matrix(new DP_entry[seq_L * parameters.m_L]);
	std::unique_ptr<T[]> right_clip_matrix(new T[seq_L - 1]);

	// initialization
	if (seq_L > 1) // if sequence is a single base it cannot be clipped
	{
		left_clip_matrix[0] = parameters.m_uniform_base_e[sequence[0]] + log_base(parameters.m_into_left_clip.from_begin);
		std::fill(right_clip_matrix.get(), right_clip_matrix.get() + (seq_L - 1), 0);
	}
	for (seq_index_type j = 1; j < seq_L - 1; ++j)
	{
		left_clip_matrix[j] = left_clip_matrix[j - 1] + parameters.m_uniform_base_e[sequence[j]] + log_base(parameters.m_into_left_clip.from_left_clip);
	}

	// first column
	for (genome_index_type i = 0; i < parameters.m_L; ++i)
	{
		DP_matrix[seq_L * (i) + (0)].ln_match = parameters.m_E[i][sequence[0]] + log_base(parameters.m_trans_matrix[i].into_match.from_begin);
		DP_matrix[seq_L * (i) + (0)].ln_insertion = std::numeric_limits<T>::lowest();
		DP_matrix[seq_L * (i) + (0)].ln_deletion = std::numeric_limits<T>::lowest();
	}

	// first row
	for (seq_index_type j = 1; j < seq_L; ++j)
	{
		DP_matrix[seq_L * (0) + (j)].ln_match = parameters.m_E[0][sequence[j]] + left_clip_matrix[j - 1] + log_base(parameters.m_trans_matrix[0].into_match.from_left_clip);
		DP_matrix[seq_L * (0) + (j)].ln_insertion = std::numeric_limits<T>::lowest();
		DP_matrix[seq_L * (0) + (j)].ln_deletion = std::numeric_limits<T>::lowest();

		right_clip_matrix[j - 1] += parameters.m_into_right_clip.from_match * exp_base(DP_matrix[seq_L * (0) + (j - 1)].ln_match);
	}

	// Calculate full DP matrix
	for (genome_index_type i = 1; i < parameters.m_L; ++i)
	{
		for (seq_index_type j = 1; j < seq_L; ++j)
		{
			// match matrix
			DP_matrix[seq_L * (i) + (j)].ln_match = parameters.m_E[i][sequence[j]] + log_base(parameters.m_trans_matrix[i].into_match.from_match * exp_base(DP_matrix[seq_L * (i - 1) + (j - 1)].ln_match) + parameters.m_trans_matrix[i].into_match.from_insertion * exp_base(DP_matrix[seq_L * (i - 1) + (j - 1)].ln_insertion) + parameters.m_trans_matrix[i].into_match.from_deletion * exp_base(DP_matrix[seq_L * (i - 1) + (j - 1)].ln_deletion) + parameters.m_trans_matrix[i].into_match.from_left_clip * exp_base(left_clip_matrix[j - 1]));

			// insertion matrix
			DP_matrix[seq_L * (i) + (j)].ln_insertion = parameters.m_uniform_base_e[sequence[j]] + log_base(parameters.m_trans_matrix[i].into_insertion.from_match * exp_base(DP_matrix[seq_L * (i) + (j - 1)].ln_match) + parameters.m_trans_matrix[i].into_insertion.from_insertion * exp_base(DP_matrix[seq_L * (i) + (j - 1)].ln_insertion));

			// deletion matrix
			DP_matrix[seq_L * (i) + (j)].ln_deletion = +log_base(parameters.m_trans_matrix[i].into_deletion.from_match * exp_base(DP_matrix[seq_L * (i - 1) + (j)].ln_match) + parameters.m_trans_matrix[i].into_deletion.from_deletion * exp_base(DP_matrix[seq_L * (i - 1) + (j)].ln_deletion));

			// right clip matrix
			right_clip_matrix[j - 1] += parameters.m_into_right_clip.from_match * exp_base(DP_matrix[seq_L * (i) + (j - 1)].ln_match);
		}
	}

	// right clip states
	if (seq_L > 1) // if sequence is a single base it cannot be clipped
	{
		right_clip_matrix[0] = parameters.m_uniform_base_e[sequence[0]] + log_base(right_clip_matrix[0]);
	}
	for (seq_index_type j = 1; j < seq_L - 1; ++j)
	{
		right_clip_matrix[j] = parameters.m_uniform_base_e[sequence[j]] + log_base(
																			  right_clip_matrix[j] + parameters.m_into_right_clip.from_right_clip * exp_base(right_clip_matrix[j - 1]));
	}

	// termination
	T result = parameters.m_into_end.from_last_match * exp_base(DP_matrix[seq_L * (parameters.m_L - 1) + (seq_L - 1)].ln_match);
	for (genome_index_type i = 0; i < parameters.m_L - 1; ++i)
	{
		result += parameters.m_into_end.from_match * exp_base(DP_matrix[seq_L * (i) + (seq_L - 1)].ln_match);
	}
	result += (seq_L > 1 ? parameters.m_into_end.from_right_clip * exp_base(right_clip_matrix[(seq_L - 1) - 1]) : 0);
	return log_base(result);
}

template <typename T>
T hmmalign<T>::viterbi(
	const reference_genome<T>& parameters,
	const boost::string_ref& sequence,
	const uint32_t ref_start,
	const uint32_t ref_end,
	std::vector<minimal_alignment>& alignments) noexcept
{
	static_assert(std::is_integral<T>::value, "T needs to be an integral type!\n");
	using genome_index_type = typename std::vector<typename reference_genome<T>::template trans_matrix<T>>::size_type;
	using seq_index_type = std::string::size_type;

	const uint32_t reference_end = std::min<uint32_t>(ref_end, parameters.m_L);

	if (ref_start >= reference_end)
	{
		std::cerr << "ERROR: Start '" << ref_start << "' has to be strictly smaller than end '" << reference_end << "'" << std::endl;
		exit(EXIT_FAILURE);
	}

	const genome_index_type region_L = reference_end - ref_start;
	const seq_index_type seq_L = sequence.length();

	/////////////////
	// DP matrices //
	/////////////////
	// left struct
	struct left_clip_struct
	{
		T score;

		bool from_begin : 1;
		bool from_left_clip : 1;
	};
	std::unique_ptr<left_clip_struct[]> left_clip_matrix(new left_clip_struct[seq_L - 1]);

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
	std::unique_ptr<DP_entry_struct[]> DP_matrix(new DP_entry_struct[seq_L * region_L]);

	// right clips
	struct right_clip_struct
	{
		T score;

		bool from_right_clip : 1;
		boost::dynamic_bitset<> from_match;
	};
	std::unique_ptr<right_clip_struct[]> right_clip_matrix(new right_clip_struct[seq_L - 1]);

	// terminal
	struct end_struct
	{
		T score;

		bool from_right_clip : 1;
		boost::dynamic_bitset<> from_match;
	} end;

	//////////////////////////
	// Phase 1: DP matrices //
	//////////////////////////
	//
	// initialization
	//
	//////////////////////
	// left clip matrix //
	//////////////////////
	if (seq_L > 1) // if sequence is a single base it cannot be clipped
	{
		left_clip_matrix[0].score = parameters.m_uniform_base_e[sequence[0]] + parameters.m_into_left_clip.from_begin;
		left_clip_matrix[0].from_begin = true;
		left_clip_matrix[0].from_left_clip = false;
	}
	for (seq_index_type j = 1; j < seq_L - 1; ++j)
	{
		left_clip_matrix[j].score = parameters.m_uniform_base_e[sequence[j]] + parameters.m_into_left_clip.from_left_clip + left_clip_matrix[j - 1].score;
		left_clip_matrix[j].from_begin = false;
		left_clip_matrix[j].from_left_clip = true;
	}

	// first column
	for (genome_index_type i = 0, genome_offset_i = ref_start; i < region_L; ++i, ++genome_offset_i)
	{
		//////////////////
		// match matrix //
		//////////////////
		DP_matrix[seq_L * (i) + (0)].match.score = parameters.m_E[ref_start + i][sequence[0]] + parameters.m_trans_matrix[ref_start + i].into_match.from_begin;
		DP_matrix[seq_L * (i) + (0)].match.from_begin = true;
		DP_matrix[seq_L * (i) + (0)].match.from_left_clip = false;
		DP_matrix[seq_L * (i) + (0)].match.from_match = false;
		DP_matrix[seq_L * (i) + (0)].match.from_insertion = false;
		DP_matrix[seq_L * (i) + (0)].match.from_deletion = false;

		//////////////////////
		// insertion matrix //
		//////////////////////
		DP_matrix[seq_L * (i) + (0)].insertion.score = lower_limit;
		DP_matrix[seq_L * (i) + (0)].insertion.from_match = false;
		DP_matrix[seq_L * (i) + (0)].insertion.from_insertion = false;

		/////////////////////
		// deletion matrix //
		/////////////////////
		DP_matrix[seq_L * (i) + (0)].deletion.score = lower_limit;
		DP_matrix[seq_L * (i) + (0)].deletion.from_match = false;
		DP_matrix[seq_L * (i) + (0)].deletion.from_deletion = false;
	}

	// first row
	for (seq_index_type j = 1; j < seq_L; ++j)
	{
		//////////////////
		// match matrix //
		//////////////////
		DP_matrix[seq_L * (0) + (j)].match.score = parameters.m_E[ref_start][sequence[j]] + parameters.m_trans_matrix[ref_start].into_match.from_left_clip + left_clip_matrix[j - 1].score;
		DP_matrix[seq_L * (0) + (j)].match.from_begin = false;
		DP_matrix[seq_L * (0) + (j)].match.from_left_clip = true;
		DP_matrix[seq_L * (0) + (j)].match.from_match = false;
		DP_matrix[seq_L * (0) + (j)].match.from_insertion = false;
		DP_matrix[seq_L * (0) + (j)].match.from_deletion = false;

		//////////////////////
		// insertion matrix //
		//////////////////////
		DP_matrix[seq_L * (0) + (j)].insertion.score = lower_limit;
		DP_matrix[seq_L * (0) + (j)].insertion.from_match = false;
		DP_matrix[seq_L * (0) + (j)].insertion.from_insertion = false;

		/////////////////////
		// deletion matrix //
		/////////////////////
		DP_matrix[seq_L * (0) + (j)].deletion.score = lower_limit;
		DP_matrix[seq_L * (0) + (j)].deletion.from_match = false;
		DP_matrix[seq_L * (0) + (j)].deletion.from_deletion = false;

		///////////////////////
		// right clip matrix //
		///////////////////////
		right_clip_matrix[j - 1].score = parameters.m_into_right_clip.from_match + DP_matrix[seq_L * (0) + (j - 1)].match.score;
		right_clip_matrix[j - 1].from_match.resize(region_L);
		right_clip_matrix[j - 1].from_match.reset();
	}

	// Calculate full DP matrix
	T temp_tr, temp_max;
	for (genome_index_type i = 1; i < region_L; ++i)
	{
		for (seq_index_type j = 1; j < seq_L; ++j)
		{
			//////////////////
			// match matrix //
			//////////////////
			DP_matrix[seq_L * (i) + (j)].match.score = parameters.m_E[ref_start + i][sequence[j]];
			DP_matrix[seq_L * (i) + (j)].match.from_begin = false;

			// match ->
			temp_max = parameters.m_trans_matrix[ref_start + i].into_match.from_match + DP_matrix[seq_L * (i - 1) + (j - 1)].match.score;
			DP_matrix[seq_L * (i) + (j)].match.from_match = true;

			// left_clip ->
			temp_tr = parameters.m_trans_matrix[ref_start + i].into_match.from_left_clip + left_clip_matrix[j - 1].score;
			if (temp_tr >= temp_max)
			{
				DP_matrix[seq_L * (i) + (j)].match.from_left_clip = true;

				if (temp_tr > temp_max)
				{
					temp_max = temp_tr;
					DP_matrix[seq_L * (i) + (j)].match.from_match = false;
				}
			}
			else
			{
				DP_matrix[seq_L * (i) + (j)].match.from_left_clip = false;
			}

			// insertion ->
			temp_tr = parameters.m_trans_matrix[ref_start + i].into_match.from_insertion + DP_matrix[seq_L * (i - 1) + (j - 1)].insertion.score;
			if (temp_tr >= temp_max)
			{
				DP_matrix[seq_L * (i) + (j)].match.from_insertion = true;

				if (temp_tr > temp_max)
				{
					temp_max = temp_tr;
					DP_matrix[seq_L * (i) + (j)].match.from_match = false;
					DP_matrix[seq_L * (i) + (j)].match.from_left_clip = false;
				}
			}
			else
			{
				DP_matrix[seq_L * (i) + (j)].match.from_insertion = false;
			}

			// deletion ->
			temp_tr = parameters.m_trans_matrix[ref_start + i].into_match.from_deletion + DP_matrix[seq_L * (i - 1) + (j - 1)].deletion.score;
			if (temp_tr >= temp_max)
			{
				DP_matrix[seq_L * (i) + (j)].match.from_deletion = true;

				if (temp_tr > temp_max)
				{
					temp_max = temp_tr;
					DP_matrix[seq_L * (i) + (j)].match.from_match = false;
					DP_matrix[seq_L * (i) + (j)].match.from_left_clip = false;
					DP_matrix[seq_L * (i) + (j)].match.from_insertion = false;
				}
			}
			else
			{
				DP_matrix[seq_L * (i) + (j)].match.from_deletion = false;
			}
			DP_matrix[seq_L * (i) + (j)].match.score += temp_max;

			//////////////////////
			// insertion matrix //
			//////////////////////
			DP_matrix[seq_L * (i) + (j)].insertion.score = parameters.m_uniform_base_e[sequence[j]];

			// match ->
			temp_max = parameters.m_trans_matrix[ref_start + i].into_insertion.from_match + DP_matrix[seq_L * (i) + (j - 1)].match.score;
			DP_matrix[seq_L * (i) + (j)].insertion.from_match = true;

			// insertion ->
			temp_tr = parameters.m_trans_matrix[ref_start + i].into_insertion.from_insertion + DP_matrix[seq_L * (i) + (j - 1)].insertion.score;
			if (temp_tr >= temp_max)
			{
				DP_matrix[seq_L * (i) + (j)].insertion.from_insertion = true;

				if (temp_tr > temp_max)
				{
					temp_max = temp_tr;
					DP_matrix[seq_L * (i) + (j)].insertion.from_match = false;
				}
			}
			else
			{
				DP_matrix[seq_L * (i) + (j)].insertion.from_insertion = false;
			}
			DP_matrix[seq_L * (i) + (j)].insertion.score += temp_max;

			/////////////////////
			// deletion matrix //
			/////////////////////
			DP_matrix[seq_L * (i) + (j)].deletion.score = 0;

			// match ->
			temp_max = parameters.m_trans_matrix[ref_start + i].into_deletion.from_match + DP_matrix[seq_L * (i - 1) + (j)].match.score;
			DP_matrix[seq_L * (i) + (j)].deletion.from_match = true;

			// deletion ->
			temp_tr = parameters.m_trans_matrix[ref_start + i].into_deletion.from_deletion + DP_matrix[seq_L * (i - 1) + (j)].deletion.score;
			if (temp_tr >= temp_max)
			{
				DP_matrix[seq_L * (i) + (j)].deletion.from_deletion = true;

				if (temp_tr > temp_max)
				{
					temp_max = temp_tr;
					DP_matrix[seq_L * (i) + (j)].deletion.from_match = false;
				}
			}
			else
			{
				DP_matrix[seq_L * (i) + (j)].deletion.from_deletion = false;
			}
			DP_matrix[seq_L * (i) + (j)].deletion.score += temp_max;

			///////////////////////
			// right clip matrix //
			///////////////////////
			temp_tr = parameters.m_into_right_clip.from_match + DP_matrix[seq_L * (i) + (j - 1)].match.score;
			if (temp_tr >= right_clip_matrix[j - 1].score)
			{
				if (temp_tr > right_clip_matrix[j - 1].score)
				{
					right_clip_matrix[j - 1].from_match.reset();
					right_clip_matrix[j - 1].score = temp_tr;
				}

				right_clip_matrix[j - 1].from_match.set(i);
			}
		}
	}

	///////////////////////
	// right clip matrix //
	///////////////////////
	for (seq_index_type j = 0; j < seq_L - 1; ++j)
	{
		temp_tr = (j > 0 ? parameters.m_into_right_clip.from_right_clip + right_clip_matrix[j - 1].score : lower_limit);
		if (temp_tr >= right_clip_matrix[j].score)
		{
			if (temp_tr > right_clip_matrix[j].score)
			{
				right_clip_matrix[j].from_match.reset();
				right_clip_matrix[j].score = temp_tr;
			}

			right_clip_matrix[j].from_right_clip = true;
		}
		else
		{
			right_clip_matrix[j].from_right_clip = false;
		}
		right_clip_matrix[j].score += parameters.m_uniform_base_e[sequence[j]];
	}

	/////////
	// end //
	/////////
	end.from_match.resize(region_L);
	end.from_match.reset();

	end.from_right_clip = true;
	end.score = (seq_L > 1 ? parameters.m_into_end.from_right_clip + right_clip_matrix[(seq_L - 1) - 1].score : lower_limit);
	for (genome_index_type i = 0; i < region_L; ++i)
	{
		temp_tr = (ref_start + i < parameters.m_L - 1 ? parameters.m_into_end.from_match : parameters.m_into_end.from_last_match) + DP_matrix[seq_L * (i) + (seq_L - 1)].match.score;

		if (temp_tr >= end.score)
		{
			if (temp_tr > end.score)
			{
				end.from_right_clip = false;

				end.from_match.reset();
				end.score = temp_tr;
			}

			end.from_match.set(i);
		}
	}

	////////////////////////
	// Phase 2: backtrack //
	////////////////////////
	alignments.clear();

	std::vector<char> cur_path;
	cur_path.reserve(1.5 * seq_L);

	struct node
	{
		char s;
		genome_index_type i;
		seq_index_type j;
	};

	std::stack<node> S;
	S.push({ 'E', 0, 0 });

	node v;
	std::size_t i;

	while (!S.empty())
	{
		v = S.top();
		S.pop();

		if (v.s == '\0')
		{
			// '\0' is a sentinel value, used
			// for when going back in the graph
			cur_path.pop_back();
		}
		else
		{
			if (v.s == 'B')
			{
				// found an alignment
				alignments.emplace_back(cur_path, v.i + ref_start);
			}
			else
			{
				// haven't found an alignment yet
				cur_path.push_back(v.s);

				S.push({ '\0', 0, 0 });

				switch (v.s)
				{
					// left clip
					case 'L':
						if (left_clip_matrix[v.j].from_begin)
						{
							S.push({ 'B', v.i, 0 });
						}
						if (left_clip_matrix[v.j].from_left_clip)
						{
							S.push({ 'L', v.i, v.j - 1 });
						}
						break;

					// match
					case 'M':
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_begin)
						{
							S.push({ 'B', v.i, 0 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_left_clip)
						{
							S.push({ 'L', v.i, v.j - 1 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_match)
						{
							S.push({ 'M', v.i - 1, v.j - 1 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_insertion)
						{
							S.push({ 'I', v.i - 1, v.j - 1 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_deletion)
						{
							S.push({ 'D', v.i - 1, v.j - 1 });
						}
						break;

					// insertion
					case 'I':
						if (DP_matrix[seq_L * (v.i) + (v.j)].insertion.from_match)
						{
							S.push({ 'M', v.i, v.j - 1 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].insertion.from_insertion)
						{
							S.push({ 'I', v.i, v.j - 1 });
						}
						break;

					// deletion
					case 'D':
						if (DP_matrix[seq_L * (v.i) + (v.j)].deletion.from_match)
						{
							S.push({ 'M', v.i - 1, v.j });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].deletion.from_deletion)
						{
							S.push({ 'D', v.i - 1, v.j });
						}
						break;

					// right clip
					case 'R':
						if (right_clip_matrix[v.j].from_right_clip)
						{
							S.push({ 'R', 0, v.j - 1 });
						}

						i = right_clip_matrix[v.j].from_match.find_first();
						while (i != boost::dynamic_bitset<>::npos)
						{
							S.push({ 'M', i, v.j });
							i = right_clip_matrix[v.j].from_match.find_next(i);
						}
						break;

					// end
					case 'E':
						if (end.from_right_clip)
						{
							S.push({ 'R', 0, seq_L - 2 });
						}

						i = end.from_match.find_first();
						while (i != boost::dynamic_bitset<>::npos)
						{
							S.push({ 'M', i, seq_L - 1 });
							i = end.from_match.find_next(i);
						}
						break;

					default:
						std::cerr << "ERROR: Unexpected state " << v.s << " while backtracking." << std::endl;
						exit(EXIT_FAILURE);
						break;
				}
			}
		}
	}

	return end.score;
}
} // unnamed namespace

#endif /* NGSHMMALIGN_HMMALIGN_IMPL_HPP */
