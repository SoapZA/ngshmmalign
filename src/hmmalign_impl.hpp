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

#include <boost/dynamic_bitset.hpp>

#include "dna_array.hpp"
#include "parameter_pack.hpp"
#include "sam.hpp"

template <typename T>
hmmalign<T>::hmmalign(const parameter_pack<T>& p_)
	: m_parameters(p_)
{
}

template <typename T>
hmmalign<T>::hmmalign(parameter_pack<T>&& p_)
	: m_parameters(std::move(p_))
{
}

template <typename T>
template <typename V>
hmmalign<T>::hmmalign(
	const std::vector<dna_array<V, 5>>& allel_freq_,
	const std::vector<V>& vec_M_to_D_p_,
	const std::vector<V>& vec_D_to_D_p_,
	const background_rates<V>& error_rates)
{
	static_assert(std::is_floating_point<V>::value, "V needs to be a floating point type!\n");
	m_parameters.set_parameters(
		allel_freq_,
		vec_M_to_D_p_,
		vec_D_to_D_p_,
		error_rates);
}

template <typename T>
template <typename V>
void hmmalign<T>::set_parameters(
	const std::vector<dna_array<V, 5>>& allel_freq_,
	const std::vector<V>& vec_M_to_D_p_,
	const std::vector<V>& vec_D_to_D_p_,
	const background_rates<V>& error_rates)
{
	static_assert(std::is_floating_point<V>::value, "V needs to be a floating point type!\n");
	m_parameters.set_parameters(
		allel_freq_,
		vec_M_to_D_p_,
		vec_D_to_D_p_,
		error_rates);
}

template <typename T>
int hmmalign<T>::get_length_profile() const
{
	return m_parameters.m_L;
}

template <typename T>
T hmmalign<T>::logLik(const std::string& sequence) const
{
	static_assert(std::is_floating_point<T>::value, "T needs to be a floating point type!\n");
	const auto seq_L = sequence.length();

	// DP matrices
	struct DP_entry
	{
		T ln_match;
		T ln_insertion;
		T ln_deletion;
	};
	std::unique_ptr<T[]> left_clip_matrix(new T[seq_L - 1]);
	std::unique_ptr<DP_entry[]> DP_matrix(new DP_entry[seq_L * m_parameters.m_L]);
	std::unique_ptr<T[]> right_clip_matrix(new T[seq_L - 1]);

	// initialization
	if (seq_L > 1) // if sequence is a single base it cannot be clipped
	{
		left_clip_matrix[0] = m_parameters.m_uniform_base_e[sequence[0]] + log_base(m_parameters.m_into_left_clip.from_begin);
		std::fill(right_clip_matrix.get(), right_clip_matrix.get() + (seq_L - 1), 0);
	}
	for (std::string::size_type j = 1; j < seq_L - 1; ++j)
	{
		left_clip_matrix[j] = left_clip_matrix[j - 1] + m_parameters.m_uniform_base_e[sequence[j]] + log_base(m_parameters.m_into_left_clip.from_left_clip);
	}

	// first column
	for (typename std::vector<typename parameter_pack<T>::template trans_matrix<T>>::size_type i = 0; i < m_parameters.m_L; ++i)
	{
		DP_matrix[seq_L * (i) + (0)].ln_match = m_parameters.m_E[i][sequence[0]] + log_base(m_parameters.m_trans_matrix[i].into_match.from_begin);
		DP_matrix[seq_L * (i) + (0)].ln_insertion = std::numeric_limits<T>::lowest();
		DP_matrix[seq_L * (i) + (0)].ln_deletion = std::numeric_limits<T>::lowest();
	}

	// first row
	for (std::string::size_type j = 1; j < seq_L; ++j)
	{
		DP_matrix[seq_L * (0) + (j)].ln_match = m_parameters.m_E[0][sequence[j]] + left_clip_matrix[j - 1] + log_base(m_parameters.m_trans_matrix[0].into_match.from_left_clip);
		DP_matrix[seq_L * (0) + (j)].ln_insertion = std::numeric_limits<T>::lowest();
		DP_matrix[seq_L * (0) + (j)].ln_deletion = std::numeric_limits<T>::lowest();

		right_clip_matrix[j - 1] += m_parameters.m_into_right_clip.from_match * exp_base(DP_matrix[seq_L * (0) + (j - 1)].ln_match);
	}

	// Calculate full DP matrix
	for (typename std::vector<typename parameter_pack<T>::template trans_matrix<T>>::size_type i = 1; i < m_parameters.m_L; ++i)
	{
		for (std::string::size_type j = 1; j < seq_L; ++j)
		{
			// match matrix
			DP_matrix[seq_L * (i) + (j)].ln_match = m_parameters.m_E[i][sequence[j]] + log_base(m_parameters.m_trans_matrix[i].into_match.from_match * exp_base(DP_matrix[seq_L * (i - 1) + (j - 1)].ln_match) + m_parameters.m_trans_matrix[i].into_match.from_insertion * exp_base(DP_matrix[seq_L * (i - 1) + (j - 1)].ln_insertion) + m_parameters.m_trans_matrix[i].into_match.from_deletion * exp_base(DP_matrix[seq_L * (i - 1) + (j - 1)].ln_deletion) + m_parameters.m_trans_matrix[i].into_match.from_left_clip * exp_base(left_clip_matrix[j - 1]));

			// insertion matrix
			DP_matrix[seq_L * (i) + (j)].ln_insertion = m_parameters.m_uniform_base_e[sequence[j]] + log_base(m_parameters.m_trans_matrix[i].into_insertion.from_match * exp_base(DP_matrix[seq_L * (i) + (j - 1)].ln_match) + m_parameters.m_trans_matrix[i].into_insertion.from_insertion * exp_base(DP_matrix[seq_L * (i) + (j - 1)].ln_insertion));

			// deletion matrix
			DP_matrix[seq_L * (i) + (j)].ln_deletion = +log_base(m_parameters.m_trans_matrix[i].into_deletion.from_match * exp_base(DP_matrix[seq_L * (i - 1) + (j)].ln_match) + m_parameters.m_trans_matrix[i].into_deletion.from_deletion * exp_base(DP_matrix[seq_L * (i - 1) + (j)].ln_deletion));

			// right clip matrix
			right_clip_matrix[j - 1] += m_parameters.m_into_right_clip.from_match * exp_base(DP_matrix[seq_L * (i) + (j - 1)].ln_match);
		}
	}

	// right clip states
	if (seq_L > 1) // if sequence is a single base it cannot be clipped
	{
		right_clip_matrix[0] = m_parameters.m_uniform_base_e[sequence[0]] + log_base(right_clip_matrix[0]);
	}
	for (std::string::size_type j = 1; j < seq_L - 1; ++j)
	{
		right_clip_matrix[j] = m_parameters.m_uniform_base_e[sequence[j]] + log_base(
																				right_clip_matrix[j] + m_parameters.m_into_right_clip.from_right_clip * exp_base(right_clip_matrix[j - 1]));
	}

	// termination
	T result = m_parameters.m_into_end.from_last_match * exp_base(DP_matrix[seq_L * (m_parameters.m_L - 1) + (seq_L - 1)].ln_match);
	for (typename std::vector<typename parameter_pack<T>::template trans_matrix<T>>::size_type i = 0; i < m_parameters.m_L - 1; ++i)
	{
		result += m_parameters.m_into_end.from_match * exp_base(DP_matrix[seq_L * (i) + (seq_L - 1)].ln_match);
	}
	result += (seq_L > 1 ? m_parameters.m_into_end.from_right_clip * exp_base(right_clip_matrix[(seq_L - 1) - 1]) : 0);
	return log_base(result);
}

template <typename T>
sam_entry hmmalign<T>::viterbi(const std::string& sequence, std::default_random_engine& generator, char clip_char) const
{
	static_assert(std::is_integral<T>::value, "T needs to be an integral type!\n");
	const auto seq_L = sequence.length();
	const T lower_limit = -1000000;

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
	std::unique_ptr<DP_entry_struct[]> DP_matrix(new DP_entry_struct[seq_L * m_parameters.m_L]);

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
		left_clip_matrix[0].score = m_parameters.m_uniform_base_e[sequence[0]] + m_parameters.m_into_left_clip.from_begin;
		left_clip_matrix[0].from_begin = true;
		left_clip_matrix[0].from_left_clip = false;
	}
	for (std::string::size_type j = 1; j < seq_L - 1; ++j)
	{
		left_clip_matrix[j].score = m_parameters.m_uniform_base_e[sequence[j]] + m_parameters.m_into_left_clip.from_left_clip + left_clip_matrix[j - 1].score;
		left_clip_matrix[j].from_begin = false;
		left_clip_matrix[j].from_left_clip = true;
	}

	// first column
	for (typename std::vector<typename parameter_pack<T>::template trans_matrix<T>>::size_type i = 0; i < m_parameters.m_L; ++i)
	{
		//////////////////
		// match matrix //
		//////////////////
		DP_matrix[seq_L * (i) + (0)].match.score = m_parameters.m_E[i][sequence[0]] + m_parameters.m_trans_matrix[i].into_match.from_begin;
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
	for (std::string::size_type j = 1; j < seq_L; ++j)
	{
		//////////////////
		// match matrix //
		//////////////////
		DP_matrix[seq_L * (0) + (j)].match.score = m_parameters.m_E[0][sequence[j]] + m_parameters.m_trans_matrix[0].into_match.from_left_clip + left_clip_matrix[j - 1].score;
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
		right_clip_matrix[j - 1].score = m_parameters.m_into_right_clip.from_match + DP_matrix[seq_L * (0) + (j - 1)].match.score;
		right_clip_matrix[j - 1].from_match.resize(m_parameters.m_L);
		right_clip_matrix[j - 1].from_match.reset();
	}

	// Calculate full DP matrix
	T temp_tr, temp_max;
	for (typename std::vector<typename parameter_pack<T>::template trans_matrix<T>>::size_type i = 1; i < m_parameters.m_L; ++i)
	{
		for (std::string::size_type j = 1; j < seq_L; ++j)
		{
			//////////////////
			// match matrix //
			//////////////////
			DP_matrix[seq_L * (i) + (j)].match.score = m_parameters.m_E[i][sequence[j]];
			DP_matrix[seq_L * (i) + (j)].match.from_begin = false;

			// match ->
			temp_max = m_parameters.m_trans_matrix[i].into_match.from_match + DP_matrix[seq_L * (i - 1) + (j - 1)].match.score;
			DP_matrix[seq_L * (i) + (j)].match.from_match = true;

			// left_clip ->
			temp_tr = m_parameters.m_trans_matrix[i].into_match.from_left_clip + left_clip_matrix[j - 1].score;
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
			temp_tr = m_parameters.m_trans_matrix[i].into_match.from_insertion + DP_matrix[seq_L * (i - 1) + (j - 1)].insertion.score;
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
			temp_tr = m_parameters.m_trans_matrix[i].into_match.from_deletion + DP_matrix[seq_L * (i - 1) + (j - 1)].deletion.score;
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
			DP_matrix[seq_L * (i) + (j)].insertion.score = m_parameters.m_uniform_base_e[sequence[j]];

			// match ->
			temp_max = m_parameters.m_trans_matrix[i].into_insertion.from_match + DP_matrix[seq_L * (i) + (j - 1)].match.score;
			DP_matrix[seq_L * (i) + (j)].insertion.from_match = true;

			// insertion ->
			temp_tr = m_parameters.m_trans_matrix[i].into_insertion.from_insertion + DP_matrix[seq_L * (i) + (j - 1)].insertion.score;
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
			temp_max = m_parameters.m_trans_matrix[i].into_deletion.from_match + DP_matrix[seq_L * (i - 1) + (j)].match.score;
			DP_matrix[seq_L * (i) + (j)].deletion.from_match = true;

			// deletion ->
			temp_tr = m_parameters.m_trans_matrix[i].into_deletion.from_deletion + DP_matrix[seq_L * (i - 1) + (j)].deletion.score;
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
			temp_tr = m_parameters.m_into_right_clip.from_match + DP_matrix[seq_L * (i) + (j - 1)].match.score;
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
	for (std::string::size_type j = 0; j < seq_L - 1; ++j)
	{
		temp_tr = (j > 0 ? m_parameters.m_into_right_clip.from_right_clip + right_clip_matrix[j - 1].score : lower_limit);
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
		right_clip_matrix[j].score += m_parameters.m_uniform_base_e[sequence[j]];
	}

	/////////
	// end //
	/////////
	end.from_match.resize(m_parameters.m_L);
	end.from_match.reset();

	end.from_right_clip = true;
	end.score = (seq_L > 1 ? m_parameters.m_into_end.from_right_clip + right_clip_matrix[(seq_L - 1) - 1].score : lower_limit);
	for (typename std::vector<typename parameter_pack<T>::template trans_matrix<T>>::size_type i = 0; i < m_parameters.m_L; ++i)
	{
		temp_tr = (i < m_parameters.m_L - 1 ? m_parameters.m_into_end.from_match : m_parameters.m_into_end.from_last_match) + DP_matrix[seq_L * (i) + (seq_L - 1)].match.score;

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
		typename std::vector<typename parameter_pack<T>::template trans_matrix<T>>::size_type i;
		std::string::size_type j;
	};

	std::stack<node> S;
	std::vector<node> cur_path;
	node v;
	std::size_t i;

	struct ext_CIGAR
	{
		std::string CIGAR;
		uint32_t POS;
		std::size_t left_clip_length = 0;
		std::size_t right_clip_length = 0;
		std::size_t segment_length = 0;

		ext_CIGAR(const std::vector<node>& alignment, uint32_t POS_, char clip_char)
			: POS(POS_)
		{
			std::stringstream ssCIGAR;

			state op = alignment.back().s;
			std::size_t count = 0;

			auto end = alignment.crend();
			for (auto it = alignment.crbegin(); it != end; ++it)
			{
				if (it->s != op)
				{
					ssCIGAR << count;
					switch (op)
					{
						case S_L:
							ssCIGAR << clip_char;
							left_clip_length = count;
							break;
						case S_R:
							ssCIGAR << clip_char;
							right_clip_length = count;
							break;
						case M:
							ssCIGAR << 'M';
							segment_length += count;
							break;
						case I:
							ssCIGAR << 'I';
							break;
						case D:
							ssCIGAR << 'D';
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

			CIGAR = ssCIGAR.str();
		}
	};

	std::vector<ext_CIGAR> all_optimal_alignments;
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
				//std::cout << "Found path starting at " << v.i << "!\n";
				//std::cout << convert2CIGAR(cur_path) << '\n';
				all_optimal_alignments.emplace_back(cur_path, v.i, clip_char);
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

	return {
		all_optimal_alignments[rand_alignment].POS + 1,
		std::move(all_optimal_alignments[rand_alignment].CIGAR),
		end.score,
		sequence.length() - all_optimal_alignments[rand_alignment].left_clip_length - all_optimal_alignments[rand_alignment].right_clip_length,
		all_optimal_alignments[rand_alignment].segment_length,
		all_optimal_alignments[rand_alignment].left_clip_length,
		all_optimal_alignments[rand_alignment].right_clip_length
	};
}

#endif /* HMMALIGN_IMPL_HPP */
