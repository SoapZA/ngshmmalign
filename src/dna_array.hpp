#ifndef NGSHMMALIGN_DNA_ARRAY_HPP
#define NGSHMMALIGN_DNA_ARRAY_HPP

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
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

template <typename T, std::size_t N>
class dna_array
{
private:
	using array_type = std::array<T, N>;
	array_type m_array;

public:
	using size_type = typename array_type::size_type;
	using value_type = typename array_type::value_type;
	using reference = value_type&;
	using const_reference = const value_type&;

	// rule-of-5 ctors
	dna_array() = default;
	dna_array(const dna_array& other) = default;
	dna_array(dna_array&& other) = default;
	dna_array& operator=(const dna_array& other) = default;
	dna_array& operator=(dna_array&& other) = default;

	// conversion ctor
	template <typename V, std::size_t M>
	friend class dna_array;

	template <typename V>
	dna_array(const dna_array<V, N>& v) noexcept
	{
		std::transform(
			v.m_array.begin(),
			v.m_array.end(),
			this->m_array.begin(),
			[](const V& i) {
				return static_cast<T>(i);
			});
	}
	template <typename V, std::size_t M>
	dna_array(dna_array<V, M>& v) noexcept : dna_array(static_cast<const dna_array<V, M>&>(v)) {}

	// aggregate initialization
	template <typename... U>
	dna_array(U&&... u) noexcept : m_array{ { std::forward<U>(u)... } } {}

	// normal indexing
	inline const_reference operator[](size_type) const noexcept;
	inline reference operator[](size_type) noexcept;

	// DNA base indexing
	inline const_reference operator[](char) const noexcept;
	inline reference operator[](char) noexcept;

	void pretty_print(std::ostream& output) const
	{
		constexpr int col_width = 15;

		output
			<< std::left << std::setw(3) << "A:" << std::left << std::setw(col_width + 1) << (*this)['A']
			<< std::left << std::setw(3) << "C:" << std::left << std::setw(col_width + 1) << (*this)['C']
			<< std::left << std::setw(3) << "G:" << std::left << std::setw(col_width + 1) << (*this)['G']
			<< std::left << std::setw(3) << "T:" << std::left << std::setw(col_width + 1) << (*this)['T']
			<< std::left << std::setw(3) << "N:" << std::left << std::setw(col_width + 1) << (*this)['N'];
	}

	// output
	friend std::ostream& operator<<(std::ostream& output, const dna_array& dna_array_) noexcept
	{
		output << dna_array_[static_cast<size_type>(0)];
		for (size_type i = 1; i < N; ++i)
		{
			output << ' ' << dna_array_[static_cast<size_type>(i)];
		}
		return output;
	}

	// input
	friend std::istream& operator>>(std::istream& input, dna_array& dna_array_) noexcept
	{
		input >> dna_array_[static_cast<size_type>(0)];
		for (size_type i = 1; i < N; ++i)
		{
			input >> dna_array_[static_cast<size_type>(i)];
		}
		return input;
	}
};

// normal indexing
template <typename T, std::size_t N>
inline typename dna_array<T, N>::const_reference dna_array<T, N>::operator[](size_type pos) const noexcept
{
	assert((0 <= pos) && (pos < N));
	return m_array[pos];
}

template <typename T, std::size_t N>
inline typename dna_array<T, N>::reference dna_array<T, N>::operator[](size_type pos) noexcept
{
	return const_cast<reference>(static_cast<const dna_array<T, N>&>(*this)[pos]);
}

// DNA base indexing
template <typename T, std::size_t N>
inline typename dna_array<T, N>::const_reference dna_array<T, N>::operator[](char base) const noexcept
{
	static_assert((N == 2) || (N == 5), "N has to be either 2 or 5! Exiting...\n");
	switch (N)
	{
		case 2:
			switch (base)
			{
				case 'A':
				case 'C':
				case 'G':
				case 'T':
					return m_array[static_cast<size_type>(0)];
				case 'N': // ambiguous base
					return m_array[static_cast<size_type>(1)];
				default:
					throw std::invalid_argument(std::string(1, base) + std::string(" is an unknown base! Exiting...\n"));
			}
			break;

		case 5:
			switch (base)
			{
				case 'A':
					return m_array[static_cast<size_type>(0)];
				case 'C':
					return m_array[static_cast<size_type>(1)];
				case 'G':
					return m_array[static_cast<size_type>(2)];
				case 'T':
					return m_array[static_cast<size_type>(3)];
				case 'N': // ambiguous base
					return m_array[static_cast<size_type>(4)];
				default:
					throw std::invalid_argument(std::string(1, base) + std::string(" is an unknown base! Exiting...\n"));
			}
			break;
	}

	// added dead code return path here, to silence ICC
	return m_array[static_cast<size_type>(0)];
}

template <typename T, std::size_t N>
inline typename dna_array<T, N>::reference dna_array<T, N>::operator[](char base) noexcept
{
	return const_cast<reference>(static_cast<const dna_array<T, N>&>(*this)[base]);
}
} // unnamed namespace

#endif /* NGSHMMALIGN_DNA_ARRAY_HPP */