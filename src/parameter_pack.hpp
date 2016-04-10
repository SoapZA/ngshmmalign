#ifndef PARAMETER_PACK_HPP
#define PARAMETER_PACK_HPP

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

#include <cassert>
#include <algorithm>
#include <clocale>
#include <cmath>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include "utility_functions.hpp"
#include "dna_array.hpp"

template <typename T>
struct background_rates
{
	// substitution rate of sequencing process
	T error_rate;
	// M->D:
	T gap_open;
	// D->D:
	T gap_extend;
	// M->I:
	T insert_open;
	// I->I:
	T insert_extend;
	// M->End:
	T end_prob;

	// Begin->S_L:
	T left_clip_open;
	// S_L->S_L:
	T left_clip_extend;
	// M->S_R:
	T right_clip_open;
	// S_R->S_R:
	T right_clip_extend;

	// output
	friend std::ostream& operator<<(std::ostream& output, const background_rates& bg_rates) noexcept
	{
		return output << std::setprecision(2)
					  << "\tError rate:              " << bg_rates.error_rate << '\n'
					  << "\tGap open:                " << bg_rates.gap_open << '\n'
					  << "\tGap extend:              " << bg_rates.gap_extend << '\n'
					  << "\tInsert open:             " << bg_rates.insert_open << '\n'
					  << "\tInsert extend:           " << bg_rates.insert_extend << '\n'
					  << "\tEnd prob:                " << bg_rates.end_prob << '\n'
					  << "\tLeft clip open:          " << bg_rates.left_clip_open << '\n'
					  << "\tLeft clip extend:        " << bg_rates.left_clip_extend << '\n'
					  << "\tRight clip open:         " << bg_rates.right_clip_open << '\n'
					  << "\tRight clip extend:       " << bg_rates.right_clip_extend << '\n';
	}
};

template <typename T>
struct parameter_pack
{
	/* length of HMM */
	std::size_t m_L;

	/* clip transition probabilities */
	// left-clip
	struct into_left_clip_struct
	{
		T from_begin;
		T from_left_clip;
	} m_into_left_clip;

	// right-clip
	struct into_right_clip_struct
	{
		T from_match;
		T from_right_clip;
	} m_into_right_clip;

	/* HMM transition probabilities */
	// complete transition probabilities
	template <typename U>
	struct trans_matrix
	{
		// match transition probabilities
		struct into_match_struct
		{
			U from_begin;
			U from_left_clip;
			U from_match;
			U from_insertion;
			U from_deletion;
		} into_match;

		// insertion transition probabilities
		struct into_insertion_struct
		{
			U from_match;
			U from_insertion;
		} into_insertion;

		// deletion transition probabilities
		struct into_deletion_struct
		{
			U from_match;
			U from_deletion;
		} into_deletion;
	};
	std::vector<trans_matrix<T>> m_trans_matrix;

	/* terminal transition probabilities */
	struct into_end_struct
	{
		T from_right_clip;
		T from_match;
		T from_last_match;
	} m_into_end;

	/* emission probabilities */
	std::vector<dna_array<T, 5>> m_E;
	dna_array<T, 2> m_uniform_base_e;

	parameter_pack() = default;

	// conversion constructors
	template <typename V>
	parameter_pack(const parameter_pack<V>&);

	template <typename V>
	void set_parameters(
		const std::vector<dna_array<V, 5>>& allel_freq_,
		const std::vector<V>& vec_M_to_D_p_,
		const std::vector<V>& vec_D_to_D_p_,
		const background_rates<V>& error_rates);

	bool display_parameters(bool fail_on_non_summation = true);
};

template <typename T>
template <typename V>
parameter_pack<T>::parameter_pack(const parameter_pack<V>& v)
{
	T(*cast)
	(const V&) = type_caster<V, T>;

	this->m_L = v.m_L;

	/* clip transition probabilities */
	// left-clip
	this->m_into_left_clip.from_begin = cast(v.m_into_left_clip.from_begin);
	this->m_into_left_clip.from_left_clip = cast(v.m_into_left_clip.from_left_clip);

	// right-clip
	this->m_into_right_clip.from_match = cast(v.m_into_right_clip.from_match);
	this->m_into_right_clip.from_right_clip = cast(v.m_into_right_clip.from_right_clip);

	/* HMM transition probabilities */
	this->m_trans_matrix.resize(m_L);
	std::transform(
		v.m_trans_matrix.begin(),
		v.m_trans_matrix.end(),
		this->m_trans_matrix.begin(),
		[cast](const typename parameter_pack<V>::template trans_matrix<V>& i) -> typename parameter_pack<T>::template trans_matrix<T>
		{
			return {
				{ cast(i.into_match.from_begin),
					cast(i.into_match.from_left_clip),
					cast(i.into_match.from_match),
					cast(i.into_match.from_insertion),
					cast(i.into_match.from_deletion) },
				{ cast(i.into_insertion.from_match),
					cast(i.into_insertion.from_insertion) },
				{ cast(i.into_deletion.from_match),
					cast(i.into_deletion.from_deletion) }
			};
		});

	/* terminal transition probabilities */
	this->m_into_end.from_right_clip = cast(v.m_into_end.from_right_clip);
	this->m_into_end.from_match = cast(v.m_into_end.from_match);
	this->m_into_end.from_last_match = cast(v.m_into_end.from_last_match);

	/* emission probabilities */
	std::vector<dna_array<T, 5>> m_E;
	this->m_E.resize(m_L);
	std::copy(v.m_E.begin(), v.m_E.end(), this->m_E.begin());

	this->m_uniform_base_e = v.m_uniform_base_e;
}

template <typename T>
template <typename V>
void parameter_pack<T>::set_parameters(
	const std::vector<dna_array<V, 5>>& allel_freq_,
	const std::vector<V>& vec_M_to_D_p_,
	const std::vector<V>& vec_D_to_D_p_,
	const background_rates<V>& error_rates)
{
	static_assert(std::is_floating_point<V>::value, "V needs to be a floating point type!\n");

	T(*cast)
	(const V&) = type_caster<V, T>;

	/* length of HMM */
	m_L = allel_freq_.size();

	/* clip transition probabilities */
	// left-clip
	m_into_left_clip.from_begin = cast(error_rates.left_clip_open);
	m_into_left_clip.from_left_clip = cast(error_rates.left_clip_extend);

	// right-clip
	m_into_right_clip.from_match = cast(error_rates.right_clip_open);
	m_into_right_clip.from_right_clip = cast(error_rates.right_clip_extend);

	/* terminal transition probabilities */
	m_into_end.from_right_clip = cast(1.0 - error_rates.right_clip_extend);
	m_into_end.from_match = cast(error_rates.end_prob);
	m_into_end.from_last_match = cast(1.0 - error_rates.right_clip_open);

	/* HMM transition probabilities */
	if (vec_M_to_D_p_.size() != m_L - 1)
	{
		std::cerr << "vec_M_to_D_p_ (L = " << vec_M_to_D_p_.size() << ") needs to have same length as emission matrix (L = " << m_L << ")!\n";
		std::terminate();
	}
	if (vec_D_to_D_p_.size() != m_L - 1)
	{
		std::cerr << "vec_D_to_D_p_ (L = " << vec_D_to_D_p_.size() << ") needs to have same length as emission matrix (L = " << m_L << ")!\n";
		std::terminate();
	}

	std::vector<trans_matrix<V>> float_trans_matrix(m_L);
	/* Position 0 in pHMM */
	// -> M transitions
	float_trans_matrix[0].into_match = {
		(1.0 - error_rates.left_clip_open) / m_L,
		(1.0 - error_rates.left_clip_extend) / m_L,
		0.0,
		0.0,
		0.0
	};
	// -> I transitions
	float_trans_matrix[0].into_insertion = {
		error_rates.insert_open,
		error_rates.insert_extend
	};
	// -> D transitions
	float_trans_matrix[0].into_deletion = {
		0.0,
		0.0
	};

	/* Position 1 in pHMM */
	V cur_M_to_D = (vec_M_to_D_p_[0] ? vec_M_to_D_p_[0] : error_rates.gap_open);
	V sum = error_rates.insert_open + cur_M_to_D + error_rates.right_clip_open + error_rates.end_prob;
	V M_to_M = 1.0 - sum;
	assert(M_to_M > 0);
	// -> M transitions
	float_trans_matrix[1].into_match = {
		(1.0 - error_rates.left_clip_open) / m_L,
		(1.0 - error_rates.left_clip_extend) / m_L,
		M_to_M,
		1.0 - error_rates.insert_extend,
		0.0
	};
	// -> I transitions
	float_trans_matrix[1].into_insertion = {
		error_rates.insert_open,
		error_rates.insert_extend
	};
	// -> D transitions
	float_trans_matrix[1].into_deletion = {
		cur_M_to_D,
		0.0
	};

	V cur_D_to_D;
	for (typename std::vector<trans_matrix<V>>::size_type i = 2; i < m_L - 1; ++i)
	{
		cur_M_to_D = (vec_M_to_D_p_[i - 1] ? vec_M_to_D_p_[i - 1] : error_rates.gap_open);

		// -> D transitions
		cur_D_to_D = (vec_D_to_D_p_[i - 1] ? vec_D_to_D_p_[i - 1] : error_rates.gap_extend);
		cur_D_to_D -= (cur_D_to_D == 1 ? error_rates.gap_open : 0);

		float_trans_matrix[i].into_deletion = {
			cur_M_to_D,
			cur_D_to_D
		};

		sum = error_rates.insert_open + cur_M_to_D + error_rates.right_clip_open + error_rates.end_prob;
		M_to_M = 1.0 - sum;
		assert(M_to_M > 0);

		// -> M transitions
		float_trans_matrix[i].into_match = {
			(1.0 - error_rates.left_clip_open) / m_L,
			(1.0 - error_rates.left_clip_extend) / m_L,
			M_to_M,
			1.0 - error_rates.insert_extend,
			1.0 - cur_D_to_D
		};

		// -> I transitions
		float_trans_matrix[i].into_insertion = {
			error_rates.insert_open,
			error_rates.insert_extend
		};
	}

	/* Position L-1 in pHMM */
	sum = error_rates.insert_open + error_rates.right_clip_open + error_rates.end_prob;
	M_to_M = 1.0 - sum;
	assert(M_to_M > 0);
	// -> M transitions
	float_trans_matrix[m_L - 1].into_match = {
		(1.0 - error_rates.left_clip_open) / m_L,
		(1.0 - error_rates.left_clip_extend) / m_L,
		M_to_M,
		1.0 - error_rates.insert_extend,
		1.0
	};
	// -> I transitions
	float_trans_matrix[m_L - 1].into_insertion = {
		0.0,
		0.0
	};
	// -> D transitions
	float_trans_matrix[m_L - 1].into_deletion = {
		0.0,
		0.0
	};

	/* copy contents back into member transition matrix */
	m_trans_matrix.resize(m_L);
	for (typename std::vector<trans_matrix<V>>::size_type i = 0; i < m_L; ++i)
	{
		m_trans_matrix[i] = {
			{ cast(float_trans_matrix[i].into_match.from_begin),
				cast(float_trans_matrix[i].into_match.from_left_clip),
				cast(float_trans_matrix[i].into_match.from_match),
				cast(float_trans_matrix[i].into_match.from_insertion),
				cast(float_trans_matrix[i].into_match.from_deletion) },
			{ cast(float_trans_matrix[i].into_insertion.from_match),
				cast(float_trans_matrix[i].into_insertion.from_insertion) },
			{ cast(float_trans_matrix[i].into_deletion.from_match),
				cast(float_trans_matrix[i].into_deletion.from_deletion) }
		};
	}

	/* emission probabilities */
	m_uniform_base_e['A'] = static_cast<T>(log_base(0.25));
	m_uniform_base_e['N'] = static_cast<T>(log_base(1.0));

	// set position-specific emission matrix
	m_E.resize(m_L);
	V temp;
	for (typename std::vector<dna_array<V, 5>>::size_type i = 0; i < m_L; ++i)
	{
		for (char j : { 'A', 'C', 'G', 'T' })
		{
			temp = 0;
			for (char k : { 'A', 'C', 'G', 'T' })
			{
				temp += allel_freq_[i][k] * (j == k ? 1 - error_rates.error_rate : error_rates.error_rate / 3);
			}
			m_E[i][j] = static_cast<T>(log_base(temp));
		}
		m_E[i]['N'] = static_cast<T>(log_base(1.0));
	}
}

template <typename V, typename std::enable_if<std::is_floating_point<V>::value, int>::type = 0>
void inline print_sum(V v, bool fail_on_non_summation, const char* suffix = "", const char* prefix = "")
{
	std::cout
		<< prefix << std::right << std::setw(30 - strlen(prefix)) << v << suffix;
	if ((std::fabs(v - 1.0) > 1E-8) && (fail_on_non_summation))
	{
		std::terminate();
	}
}

template <typename V, typename std::enable_if<std::is_integral<V>::value, int>::type = 0>
void inline print_sum(V v, bool fail_on_non_summation, const char* suffix = "", const char* prefix = "")
{
}

template <typename T>
bool parameter_pack<T>::display_parameters(bool fail_on_non_summation)
{
	std::setlocale(LC_ALL, "en_US.UTF-8");
	std::cout
		<< std::setprecision(3) << std::fixed;
	std::cout
		<< "HMM Transition tables for <" << typeid(m_trans_matrix[0].into_match.from_left_clip).name() << ">:\n";

	const int tr_col_width = 8;
	const int value_col_width = 9;

	/* BEGIN */
	std::cout
		<< std::left << std::setw(tr_col_width) << "BEG"
		<< u8"\u2500\u252C\u2500> "
		<< std::left << std::setw(5 + 2) << "S_L" << ':'
		<< std::right << std::setw(value_col_width) << m_into_left_clip.from_begin << '\n';
	std::cout
		<< std::left << std::string(tr_col_width, ' ')
		<< u8" \u2514\u2500> "
		<< std::left << std::setw(5 + 2) << "M_j" << ':'
		<< std::right << std::setw(value_col_width) << m_trans_matrix[0].into_match.from_begin << '\n';
	print_sum<T>(m_into_left_clip.from_begin + m_L * m_trans_matrix[0].into_match.from_begin, fail_on_non_summation, "\n", "Sum:");

	/* LEFT CLIP */
	std::cout
		<< "\n\n"
		<< std::left << std::setw(tr_col_width) << "S_L"
		<< u8"\u2500\u252C\u2500> "
		<< std::left << std::setw(5 + 2) << "S_L" << ':'
		<< std::right << std::setw(value_col_width) << m_into_left_clip.from_left_clip << '\n';
	std::cout
		<< std::left << std::string(tr_col_width, ' ')
		<< u8" \u2514\u2500> "
		<< std::left << std::setw(5 + 2) << "M_j" << ':'
		<< std::right << std::setw(value_col_width) << m_trans_matrix[0].into_match.from_left_clip << '\n';
	print_sum<T>(m_into_left_clip.from_left_clip + m_L * m_trans_matrix[0].into_match.from_left_clip, fail_on_non_summation, "\n", "Sum:");

	/* MAIN STATES */
	for (typename std::vector<trans_matrix<T>>::size_type i = 0; i < m_L - 1; ++i)
	{
		// first row
		std::cout
			<< "\n\n"
			<< std::left << "M_" << std::setw(tr_col_width - 2) << i
			<< u8"\u2500\u252C\u2500> "
			<< std::left << "M_" << std::setw(5) << i + 1 << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_match.from_match;
		std::cout
			<< std::string(2 * tr_col_width, ' ')
			<< std::left << "I_" << std::setw(tr_col_width - 2) << i
			<< u8"\u2500\u252C\u2500> "
			<< std::left << "M_" << std::setw(5) << i + 1 << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_match.from_insertion;

		if (i != 0)
		{
			std::cout
				<< std::string(2 * tr_col_width, ' ')
				<< std::left << "D_" << std::setw(tr_col_width - 2) << i
				<< u8"\u2500\u252C\u2500> "
				<< std::left << "M_" << std::setw(5) << i + 1 << ':'
				<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_match.from_deletion;
		}

		// second row
		std::cout
			<< '\n'
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u251C\u2500> "
			<< std::left << "D_" << std::setw(5) << i + 1 << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_deletion.from_match;
		std::cout
			<< std::string(2 * tr_col_width, ' ')
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u2514\u2500> "
			<< std::left << "I_" << std::setw(5) << i << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i].into_insertion.from_insertion;

		if (i != 0)
		{
			std::cout
				<< std::string(2 * tr_col_width, ' ')
				<< std::left << std::string(tr_col_width, ' ')
				<< u8" \u2514\u2500> "
				<< std::left << "D_" << std::setw(5) << i + 1 << ':'
				<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_deletion.from_deletion;
		}

		// third row
		std::cout
			<< '\n'
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u251C\u2500> "
			<< std::left << "I_" << std::setw(5) << i << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i].into_insertion.from_match << '\n';

		// fourth row
		std::cout
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u251C\u2500> "
			<< std::left << std::setw(5 + 2) << "S_R" << ':'
			<< std::right << std::setw(value_col_width) << m_into_right_clip.from_match << '\n';

		// fifth row
		std::cout
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u2514\u2500> "
			<< std::left << std::setw(5 + 2) << "END" << ':'
			<< std::right << std::setw(value_col_width) << m_into_end.from_match << '\n';

		print_sum<T>(
			m_trans_matrix[i + 1].into_match.from_match + m_trans_matrix[i + 1].into_deletion.from_match + m_trans_matrix[i].into_insertion.from_match + m_into_right_clip.from_match + m_into_end.from_match,
			fail_on_non_summation, "", "Sum:");

		std::cout << std::string(2 * tr_col_width, ' ');
		print_sum<T>(
			m_trans_matrix[i + 1].into_match.from_insertion + m_trans_matrix[i].into_insertion.from_insertion,
			fail_on_non_summation);
		if (i != 0)
		{
			std::cout << std::string(2 * tr_col_width, ' ');
			print_sum<T>(
				m_trans_matrix[i + 1].into_match.from_deletion + m_trans_matrix[i + 1].into_deletion.from_deletion,
				fail_on_non_summation, "\n");
		}
	}

	/* LAST MATCH STATE */
	std::cout
		<< "\n\n"
		<< std::left << "M_" << std::setw(tr_col_width - 2) << m_L - 1
		<< u8"\u2500\u252C\u2500> "
		<< std::left << std::setw(5 + 2) << "S_R" << ':'
		<< std::right << std::setw(value_col_width) << m_into_right_clip.from_match << "\n";
	std::cout
		<< std::left << std::string(tr_col_width, ' ')
		<< u8" \u2514\u2500> "
		<< std::left << std::setw(5 + 2) << "END" << ':'
		<< std::right << std::setw(value_col_width) << m_into_end.from_last_match << '\n';
	print_sum<T>(m_into_right_clip.from_match + m_into_end.from_last_match, fail_on_non_summation, "\n", "Sum:");

	/* RIGHT CLIP */
	std::cout
		<< "\n\n"
		<< std::left << std::setw(tr_col_width) << "S_R"
		<< u8"\u2500\u252C\u2500> "
		<< std::left << std::setw(5 + 2) << "S_R" << ':'
		<< std::right << std::setw(value_col_width) << m_into_right_clip.from_right_clip << '\n';
	std::cout
		<< std::left << std::string(tr_col_width, ' ')
		<< u8" \u2514\u2500> "
		<< std::left << std::setw(5 + 2) << "END" << ':'
		<< std::right << std::setw(value_col_width) << m_into_end.from_right_clip << '\n';
	print_sum<T>(m_into_right_clip.from_right_clip + m_into_end.from_right_clip, fail_on_non_summation, "\n", "Sum:");

	int j = 0;
	std::cout << "\nEmission tables (log):\n";
	for (const auto& i : m_E)
	{
		std::cout << j++ << ":\t" << i;
	}
	std::cout << "\nEmission tables (uniform):\n\t" << m_uniform_base_e;

	return true;
}

#endif /* PARAMETER_PACK_HPP */
