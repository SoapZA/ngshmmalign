#ifndef HMMALIGN_IMPL_HPP
#define HMMALIGN_IMPL_HPP

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

#include <vector>
#include <string>
#include <sstream>
#include <list>
#include <stack>
#include <utility>
#include <cmath>
#include <random>
#include <algorithm>

#include <boost/dynamic_bitset.hpp>

#include "dna_array.hpp"
#include "sam.hpp"

namespace
{

template <typename T>
T hmmalign<T>::logLik(
	const reference_genome<T>& parameters,
	const boost::string_ref& sequence) const
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
sam_entry hmmalign<T>::viterbi(
	const reference_genome<T>& parameters,
	const boost::string_ref& sequence,
	std::default_random_engine& generator,
	const uint32_t ref_start,
	const uint32_t ref_end,
	const bool differentiate_match_state) const
{
	static_assert(std::is_integral<T>::value, "T needs to be an integral type!\n");
	using genome_index_type = typename std::vector<typename reference_genome<T>::template trans_matrix<T>>::size_type;
	using seq_index_type = std::string::size_type;

	const uint32_t reference_end = std::min<uint32_t>(ref_end, parameters.m_L);

	if (ref_start >= reference_end)
	{
		std::cerr << "ERROR: Start '" << ref_start << "' has to be strictly smaller than end '" << reference_end << "'\n";
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
	enum state
	{
		B,
		S_L,
		M,
		I,
		D,
		S_R,
		E,
		sentinel
	};

	struct node
	{
		state s;
		genome_index_type i;
		seq_index_type j;
	};

	std::stack<node> S;
	std::vector<node> cur_path;
	node v;
	std::size_t i;

	struct alignment_record
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

		CIGAR_vec CIGAR;

		uint32_t POS;
		uint32_t left_clip_length = 0;
		uint32_t right_clip_length = 0;
		uint32_t segment_length = 0;

		uint32_t edit_distance = 0;

		alignment_record(const std::vector<node>& alignment, uint32_t POS_)
			: POS(POS_)
		{
			CIGAR.reserve(5);

			state op = alignment.back().s;
			uint32_t count = 0;

			const auto it_end = alignment.crend();
			for (auto it = alignment.crbegin(); it != it_end; ++it)
			{
				if (it->s != op)
				{
					switch (op)
					{
						case S_L:
							left_clip_length = count;
							break;
						case S_R:
							right_clip_length = count;
							break;
						case M:
							CIGAR.emplace_back('M', count);
							segment_length += count;
							break;
						case I:
							CIGAR.emplace_back('I', count);
							break;
						case D:
							CIGAR.emplace_back('D', count);
							segment_length += count;
							break;
						default:
							break;
					}
					op = it->s;
					count = 1;
				}
				else
				{
					++count;
				}
			}

			//CIGAR = ssCIGAR.str();
		}
	};

	std::vector<alignment_record> all_optimal_alignments;
	S.push({ E, 0, 0 });

	while (!S.empty())
	{
		v = S.top();
		S.pop();

		if (v.s == sentinel)
		{
			// going back in the graph
			cur_path.pop_back();
		}
		else
		{
			if (v.s == B)
			{
				// found an alignment
				all_optimal_alignments.emplace_back(cur_path, v.i + ref_start);
			}
			else
			{
				// haven't found an alignment yet
				cur_path.push_back(v);

				S.push({ sentinel, 0, 0 });

				switch (v.s)
				{
					case S_L:
						if (left_clip_matrix[v.j].from_begin)
						{
							S.push({ B, v.i, 0 });
						}
						if (left_clip_matrix[v.j].from_left_clip)
						{
							S.push({ S_L, v.i, v.j - 1 });
						}
						break;

					case M:
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_begin)
						{
							S.push({ B, v.i, 0 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_left_clip)
						{
							S.push({ S_L, v.i, v.j - 1 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_match)
						{
							S.push({ M, v.i - 1, v.j - 1 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_insertion)
						{
							S.push({ I, v.i - 1, v.j - 1 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].match.from_deletion)
						{
							S.push({ D, v.i - 1, v.j - 1 });
						}
						break;

					case I:
						if (DP_matrix[seq_L * (v.i) + (v.j)].insertion.from_match)
						{
							S.push({ M, v.i, v.j - 1 });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].insertion.from_insertion)
						{
							S.push({ I, v.i, v.j - 1 });
						}
						break;

					case D:
						if (DP_matrix[seq_L * (v.i) + (v.j)].deletion.from_match)
						{
							S.push({ M, v.i - 1, v.j });
						}
						if (DP_matrix[seq_L * (v.i) + (v.j)].deletion.from_deletion)
						{
							S.push({ D, v.i - 1, v.j });
						}
						break;

					case S_R:
						if (right_clip_matrix[v.j].from_right_clip)
						{
							S.push({ S_R, 0, v.j - 1 });
						}

						i = right_clip_matrix[v.j].from_match.find_first();
						while (i != boost::dynamic_bitset<>::npos)
						{
							S.push({ M, i, v.j });
							i = right_clip_matrix[v.j].from_match.find_next(i);
						}
						break;

					case E:
						if (end.from_right_clip)
						{
							S.push({ S_R, 0, seq_L - 2 });
						}

						i = end.from_match.find_first();
						while (i != boost::dynamic_bitset<>::npos)
						{
							S.push({ M, i, seq_L - 1 });
							i = end.from_match.find_next(i);
						}
						break;

					default:
						break;
				}
			}
		}
	}

	std::uniform_int_distribution<uint32_t> range(0, all_optimal_alignments.size() - 1);
	uint32_t rand_alignment = range(generator);

	// calculate edit/Hamming distance
	struct calculate_edit_distance
	{
		void operator()(const reference_genome<T>& parameters, const boost::string_ref& seq, alignment_record& alignment, bool rewrite_cigar) const
		{
			uint32_t pos_in_read = alignment.left_clip_length;
			uint32_t pos_on_genome = alignment.POS;

			CIGAR_vec new_CIGAR_vec;
			char old_state, cur_state;
			uint32_t cur_length;

			for (const auto& op : alignment.CIGAR)
			{
				switch (op.first)
				{
					case 'M':
						if (rewrite_cigar)
						{
							old_state = (parameters.m_table_of_included_bases[pos_on_genome][seq[pos_in_read]] ? '=' : 'X');
							cur_length = 0;
						}

						for (uint32_t j = pos_in_read, i = pos_on_genome; j < pos_in_read + op.second; ++i, ++j)
						{
							cur_state = (parameters.m_table_of_included_bases[i][seq[j]] ? '=' : 'X');
							alignment.edit_distance += (cur_state == 'X');

							if (rewrite_cigar)
							{
								if (cur_state != old_state)
								{
									// switch states
									new_CIGAR_vec.emplace_back(old_state, cur_length);
									old_state = cur_state;
									cur_length = 0;
								}

								++cur_length;
							}
						}

						if (rewrite_cigar)
						{
							new_CIGAR_vec.emplace_back(cur_state, cur_length);
						}

						pos_on_genome += op.second;
						pos_in_read += op.second;
						break;

					case 'I':
						if (rewrite_cigar)
						{
							new_CIGAR_vec.push_back(op);
						}

						pos_in_read += op.second;
						break;

					case 'D':
						if (rewrite_cigar)
						{
							new_CIGAR_vec.push_back(op);
						}

						pos_on_genome += op.second;
						break;
				}
			}

			if (rewrite_cigar)
			{
				alignment.CIGAR.swap(new_CIGAR_vec);
			}
		}
	};

	// replace 'M' in CIGAR if differentiate_match_state == true
	// Match: '=' Mismatch: 'X'
	calculate_edit_distance()(parameters, sequence, all_optimal_alignments[rand_alignment], differentiate_match_state);

	return {
		parameters.m_reference_genome_name.c_str(),
		all_optimal_alignments[rand_alignment].POS,
		std::move(all_optimal_alignments[rand_alignment].CIGAR),
		end.score,
		all_optimal_alignments[rand_alignment].edit_distance,
		static_cast<int32_t>(sequence.length()) - all_optimal_alignments[rand_alignment].left_clip_length - all_optimal_alignments[rand_alignment].right_clip_length,
		all_optimal_alignments[rand_alignment].segment_length,
		all_optimal_alignments[rand_alignment].left_clip_length,
		all_optimal_alignments[rand_alignment].right_clip_length
	};
}
}

#endif /* HMMALIGN_IMPL_HPP */