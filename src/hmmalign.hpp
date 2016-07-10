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
#include "sam.hpp"
#include "reference.hpp"

namespace
{

template <typename T>
class hmmalign
{
public:
	// MARGINAL PROBABILITY of a sequence
	static T logLik(const reference_genome<T>& parameters, const boost::string_ref& sequence);

	static inline T Lik(const reference_genome<T>& parameters, const boost::string_ref& sequence)
	{
		return exp_base(logLik(parameters, sequence));
	}

	// OPTIMAL ALIGNMENT of a sequence
	static sam_entry viterbi(
		const reference_genome<T>& parameters,
		const boost::string_ref& sequence,
		std::default_random_engine& generator,
		uint32_t start,
		uint32_t end,
		bool differentiate_match_state);
};
}

#include "hmmalign_impl.hpp"

#endif /* HMMALIGN_HPP */