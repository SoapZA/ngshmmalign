#ifndef HMMALIGN_HPP
#define HMMALIGN_HPP

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

#define MAGIC_NUMBER 0xDEADBEEF

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
class hmmalign
{
public:
	// Constructors:
	hmmalign() = default;

	template <typename V>
	hmmalign(
		const std::vector<dna_array<V, 5>>& allel_freq_,
		const std::vector<V>& vec_M_to_D_p_,
		const std::vector<V>& vec_D_to_D_p_,
		const background_rates<V>& error_rates);

	hmmalign(const parameter_pack<T>&);
	hmmalign(parameter_pack<T>&&);

	hmmalign(const hmmalign& other) = default;
	hmmalign(hmmalign&& other) = default;
	hmmalign& operator=(const hmmalign& other) = default;
	hmmalign& operator=(hmmalign&& other) = default;

	template <typename V>
	friend class hmmalign;

	template <typename V>
	hmmalign(const hmmalign<V>& other_hmmalign) noexcept
	{
		m_parameters = other_hmmalign.m_parameters;
	}

	// Destructors:
	~hmmalign() = default;

	// Member functions:
	template <typename V>
	void set_parameters(
		const std::vector<dna_array<V, 5>>& allel_freq_,
		const std::vector<V>& vec_M_to_D_p_,
		const std::vector<V>& vec_D_to_D_p_,
		const background_rates<V>& error_rates);

	int get_length_profile() const;

	// MARGINAL PROBABILITY of a sequence
	T logLik(const std::string& sequence) const;

	inline T Lik(const std::string& sequence) const
	{
		return exp_base(logLik(sequence));
	}

	// OPTIMAL ALIGNMENT of a sequence
	sam_entry viterbi(const std::string& sequence, std::default_random_engine& generator, char clip_char) const;

private:
	// parameters
	parameter_pack<T> m_parameters;
};

#include "hmmalign_impl.hpp"

#endif /* HMMALIGN_HPP */
