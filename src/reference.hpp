#ifndef NGSHMMALIGN_REFERENCE_HPP
#define NGSHMMALIGN_REFERENCE_HPP

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
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <boost/utility/string_ref.hpp>

#include "dna_array.hpp"
#include "utility_functions.hpp"

extern int num_threads;

namespace
{

struct background_rates
{
	// frequency cutoff for considering a base part of the reference genome
	double low_frequency_cutoff;
	// substitution rate of sequencing process
	double error_rate;
	// M->D:
	double gap_open;
	// D->D:
	double gap_extend;
	// M->I:
	double insert_open;
	// I->I:
	double insert_extend;
	// M->End:
	double end_prob;

	// Begin->S_L:
	double left_clip_open;
	// S_L->S_L:
	double left_clip_extend;
	// M->S_R:
	double right_clip_open;
	// S_R->S_R:
	double right_clip_extend;

	// output
	friend std::ostream& operator<<(std::ostream& output, const background_rates& bg_rates) noexcept
	{
		return output << std::setprecision(2)
					  << "\tNumber of threads -t:    " << num_threads << std::endl
					  << std::endl
					  << "\tLower frequency cutoff:  " << bg_rates.low_frequency_cutoff << std::endl
					  << "\tError rate:              " << bg_rates.error_rate << std::endl
					  << "\tGap open:                " << bg_rates.gap_open << std::endl
					  << "\tGap extend:              " << bg_rates.gap_extend << std::endl
					  << "\tInsert open:             " << bg_rates.insert_open << std::endl
					  << "\tInsert extend:           " << bg_rates.insert_extend << std::endl
					  << "\tEnd prob:                " << bg_rates.end_prob << std::endl
					  << "\tLeft clip open:          " << bg_rates.left_clip_open << std::endl
					  << "\tLeft clip extend:        " << bg_rates.left_clip_extend << std::endl
					  << "\tRight clip open:         " << bg_rates.right_clip_open << std::endl
					  << "\tRight clip extend:       " << bg_rates.right_clip_extend << std::endl;
	}
};

struct reference_haplotype
{
	std::string name;
	std::string sequence;

	// we use half-open intervals
	// i.e., for the enclosed-interval [start, end)
	int32_t start;
	int32_t end;
	double count;

	reference_haplotype(const std::string& id, std::string&& seq) noexcept;
};

struct boost_string_ref_hash
{
	inline std::size_t operator()(const boost::string_ref& str_ref) const
	{
		return boost::hash_range(str_ref.begin(), str_ref.end());
	}
};

struct boost_string_ref_equal_to
{
	inline bool operator()(const boost::string_ref& str_ref1, const boost::string_ref& str_ref2) const
	{
		return str_ref1 == str_ref2;
	}
};

template <typename T>
struct reference_genome
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

	/* table of bases, without errors */
	bool m_emission_tables_initialized = false;
	std::vector<dna_array<bool, 5>> m_table_of_included_bases;

	/* majority reference sequence */
	std::string m_majority_ref;

	/* ambiguous reference sequence */
	std::string m_ambig_ref;

	uint32_t read_length_profile;

	// Ctor
	reference_genome() = default;

	template <typename V>
	reference_genome(const reference_genome<V>& v);

	// Setter
	struct msa_input
	{
	};
	struct serialized_input
	{
	};

	void set_parameters(
		std::vector<dna_array<double, 5>>&& allel_freq_,
		std::vector<double>&& vec_M_to_D_p_,
		std::vector<double>&& vec_D_to_D_p_,
		const background_rates& error_rates,
		bool ambig_bases_unequal_weight);

	void set_parameters(
		const std::string& input_msa,
		const background_rates& error_rates,
		uint32_t read_lengths,
		bool ambig_bases_unequal_weight,
		msa_input);

	void set_parameters(
		const std::string& input_hmm_archive,
		const background_rates& error_rates,
		uint32_t read_lengths,
		bool ambig_bases_unequal_weight,
		serialized_input);

	void set_parameters(
		const std::vector<reference_haplotype>& refs,
		const background_rates& error_rates,
		bool ambig_bases_unequal_weight);

	void init_emission_table(
		const std::vector<dna_array<double, 5>>& allel_freq_,
		double min_freq,
		double error_rate,
		bool ambig_bases_unequal_weight,
		bool sufficient_coverage_check);

	bool display_parameters(std::ostream& output, bool fail_on_non_summation = true) const;

	void save_to_file(std::ostream& ofs);

	// kmer-based lookup
	void create_index(uint16_t desired_kmer_length = 20);

	struct index_stat
	{
		int32_t POS;
		int32_t sd;
		uint32_t num_samples;
	};

	index_stat find_pos(const boost::string_ref& query) const;

private:
	void init_parameters(
		const background_rates& error_rates,
		bool ambig_bases_unequal_weight);

	void load_from_file(std::istream& ifs);

	// parameters from which HMM is built
	std::vector<dna_array<double, 5>> m_allel_freq;
	std::vector<double> m_vec_M_to_D_p;
	std::vector<double> m_vec_D_to_D_p;

	uint16_t m_kmer_length;
	using hash_map_type = boost::unordered_map<std::string, std::vector<uint32_t>, boost_string_ref_hash, boost_string_ref_equal_to>;
	hash_map_type m_kmer_index;
};
} // unnamed namespace

#include "index_impl.hpp"
#include "reference_impl.hpp"

#endif /* NGSHMMALIGN_REFERENCE_HPP */