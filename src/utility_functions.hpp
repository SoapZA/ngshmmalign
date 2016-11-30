#ifndef NGSHMMALIGN_UTILITY_FUNCTIONS_HPP
#define NGSHMMALIGN_UTILITY_FUNCTIONS_HPP

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
#include <map>
#include <string>
#include <type_traits>

namespace
{

static double log_exp_base = 1.10;
static short int_add = 0;
constexpr int lower_limit = -10000000;

constexpr double MAGIC_NUMBER = 0xDEADBEEF;

static const std::map<const std::string, const char> ambig_to_wobble_base = {
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
static const std::map<const char, const std::string> wobble_to_ambig_bases = {
	{ 'A', "A" },
	{ 'C', "C" },
	{ 'G', "G" },
	{ 'T', "T" },
	{ 'M', "AC" },
	{ 'R', "AG" },
	{ 'W', "AT" },
	{ 'S', "CG" },
	{ 'Y', "CT" },
	{ 'K', "GT" },
	{ 'V', "ACG" },
	{ 'H', "ACT" },
	{ 'D', "AGT" },
	{ 'B', "CGT" },
	{ 'N', "ACGT" }
};

inline char identity_char(char c) noexcept
{
	return c;
}

inline char rev_comp_char(char c) noexcept
{
	switch (c)
	{
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'N':
			return 'N';
		default:
			throw std::invalid_argument(std::string(1, c) + std::string(" is an unknown base! Exiting...\n"));
	}
}

template <typename T>
inline bool AlmostEqualRelative(T A, T B, T maxRelDiff = 1E-8)
{
	// Calculate the difference.
	T diff = std::fabs(A - B);
	A = std::fabs(A);
	B = std::fabs(B);
	// Find the largest
	T largest = std::max(A, B);

	return (diff <= largest * maxRelDiff);
}

// log function with arbitrary base
template <typename floating_point>
inline floating_point log_base(floating_point argument) noexcept
{
	static_assert(std::is_floating_point<floating_point>::value, "The log_base function requires a floating point type");
	return std::log(argument) / std::log(log_exp_base);
}

// exp function with arbitrary base
template <typename floating_point>
inline floating_point exp_base(floating_point argument) noexcept
{
	return std::exp(std::log(log_exp_base) * argument);
}

// cast from floating point -> integral
template <typename from_T, typename to_T, typename std::enable_if<std::is_floating_point<from_T>::value && std::is_integral<to_T>::value, int>::type = 0>
inline to_T type_caster(const from_T& from)
{
	if (from > 0)
	{
		return static_cast<to_T>(log_base<from_T>(from)) + int_add;
	}
	else
	{
		return lower_limit;
	}
}

// cast from integral -> floating point
template <typename from_T, typename to_T, typename std::enable_if<std::is_integral<from_T>::value && std::is_floating_point<to_T>::value, int>::type = 0>
inline to_T type_caster(const from_T& from)
{
	return static_cast<to_T>(exp_base<to_T>(from));
}

// cast from floating point -> floating point
template <typename from_T, typename to_T, typename std::enable_if<std::is_floating_point<from_T>::value && std::is_floating_point<to_T>::value, int>::type = 0>
inline to_T type_caster(const from_T& from)
{
	return static_cast<to_T>(from);
}

// cast from integral -> integral
template <typename from_T, typename to_T, typename std::enable_if<std::is_integral<from_T>::value && std::is_integral<to_T>::value, int>::type = 0>
inline to_T type_caster(const from_T& from)
{
	return static_cast<to_T>(from);
}
} // unnamed namespace

#endif /* NGSHMMALIGN_UTILITY_FUNCTIONS_HPP */
