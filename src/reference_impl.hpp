#ifndef NGSHMMALIGN_REFERENCE_IMPL_HPP
#define NGSHMMALIGN_REFERENCE_IMPL_HPP

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

#include <clocale>
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
#include "fasta.hpp"

namespace
{

reference_haplotype::reference_haplotype(const std::string& id, std::string&& seq) noexcept : sequence(std::move(seq))
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
			std::cerr << "ERROR: Count argument '" << split_vec[1] << "' is not an integral/floating point value! Aborting." << std::endl;
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

// Ctor
template <typename T>
template <typename V>
reference_genome<T>::reference_genome(const reference_genome<V>& v)
	: m_L(v.m_L),
	  m_emission_tables_initialized(v.m_emission_tables_initialized),
	  m_table_of_included_bases(v.m_table_of_included_bases),
	  m_majority_ref(v.m_majority_ref),
	  m_ambig_ref(v.m_ambig_ref),
	  read_length_profile(v.read_length_profile),
	  m_allel_freq(v.m_allel_freq),
	  m_vec_M_to_D_p(v.m_vec_M_to_D_p),
	  m_vec_D_to_D_p(v.m_vec_D_to_D_p),
	  m_kmer_length(v.m_kmer_length),
	  m_kmer_index(v.m_kmer_index)
{
	using V_ref_type = typename reference_genome<V>::template trans_matrix<V>;
	using T_ref_type = typename reference_genome<T>::template trans_matrix<T>;

	T(*cast)
	(const V&) = type_caster<V, T>;

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
		[cast](const V_ref_type& i) -> T_ref_type {
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
	this->m_E.resize(m_L);
	std::copy(v.m_E.cbegin(), v.m_E.cend(), this->m_E.begin());

	this->m_uniform_base_e = v.m_uniform_base_e;
}

// Setter
template <typename T>
void reference_genome<T>::set_parameters(
	std::vector<dna_array<double, 5>>&& allel_freq_,
	std::vector<double>&& vec_M_to_D_p_,
	std::vector<double>&& vec_D_to_D_p_,
	const background_rates& error_rates,
	bool ambig_bases_unequal_weight)
{
	m_allel_freq = std::move(allel_freq_);
	m_vec_M_to_D_p = std::move(vec_M_to_D_p_);
	m_vec_D_to_D_p = std::move(vec_D_to_D_p_);

	init_parameters(error_rates, ambig_bases_unequal_weight);
}

template <typename T>
void reference_genome<T>::init_emission_table(
	const std::vector<dna_array<double, 5>>& allel_freq_,
	const double min_freq,
	const double error_rate,
	const bool ambig_bases_unequal_weight,
	const bool sufficient_coverage_check)
{
	if (allel_freq_.size() != m_L)
	{
		std::cerr << "ERROR: New emission table 'table_base_abundance' (= " << allel_freq_.size() << ") and HMM length (= " << m_L << ") are unequal." << std::endl;
		exit(EXIT_FAILURE);
	}

	m_uniform_base_e['A'] = static_cast<T>(log_base(0.25));
	m_uniform_base_e['N'] = static_cast<T>(log_base(1.0));

	// set position-specific emission matrix
	m_allel_freq.resize(m_L);
	m_E.resize(m_L);
	m_table_of_included_bases.resize(m_L);

	// determine majority and ambiguous reference sequence
	std::string new_majority_ref;
	std::string new_ambig_ref;

	std::default_random_engine rng;

	for (typename std::vector<dna_array<double, 5>>::size_type i = 0; i < allel_freq_.size(); ++i)
	{
		dna_array<double, 5> cur_allel_table = allel_freq_[i];

		// 1. determine coverage
		const double coverage = std::accumulate(&cur_allel_table[static_cast<std::size_t>(0)], &cur_allel_table[static_cast<std::size_t>(4)], 0.0);
		if (coverage < 0)
		{
			// error out
			std::cerr << "ERROR: Coverage is negative at position " << i << " of genome." << std::endl;
			exit(EXIT_FAILURE);
		}
		else
		{
			if (coverage > 0)
			{
				// have coverage, can call a base
				// 2. determine majority + ambiguous bases, normalize frequency components
				std::string majority_bases, ambig_bases;
				double majority_freq = -1;
				double renormalization_sum = 0;
				for (char j : { 'A', 'C', 'G', 'T' })
				{
					double& allel_freq = cur_allel_table[j];
					allel_freq /= coverage;

					if (allel_freq < 0)
					{
						// error out
						std::cerr << "ERROR: Base '" << j << "' has negative coverage at position " << i << " of genome." << std::endl;
						exit(EXIT_FAILURE);
					}

					// majority base
					if (allel_freq >= majority_freq)
					{
						if (allel_freq > majority_freq)
						{
							majority_bases.clear();
							majority_freq = allel_freq;
						}

						majority_bases.push_back(j);
					}

					// ambiguous bases
					if (allel_freq > min_freq)
					{
						ambig_bases.push_back(j);
						renormalization_sum += allel_freq;
					}
					else
					{
						allel_freq = 0;
					}
				}

				// 3. call new majority + ambiguous bases, normalize frequencies for bases kept
				char new_ambig_base, new_majority_base = majority_bases[std::uniform_int_distribution<uint8_t>(0, majority_bases.length() - 1)(rng)];

				const auto it = ambig_to_wobble_base.find(ambig_bases);
				if (it != ambig_to_wobble_base.end())
				{
					new_ambig_base = it->second;
				}
				else
				{
					std::cerr << "ERROR: Could not map " << ambig_bases << " to a wobble base." << std::endl;
					exit(EXIT_FAILURE);
				}

				if ((sufficient_coverage_check) && (coverage <= std::ceil(2.0 / min_freq)))
				{
					// low coverage, use majority base instead
					cur_allel_table = { 0.0, 0.0, 0.0, 0.0, 0.0 };
					cur_allel_table[new_majority_base] = 1.0;

					new_majority_base = std::tolower(new_majority_base);
					new_ambig_base = new_majority_base;
				}
				else
				{
					// just renormalize components
					for (char j : { 'A', 'C', 'G', 'T' })
					{
						cur_allel_table[j] /= renormalization_sum;
					}
				}

				new_majority_ref.push_back(new_majority_base);
				new_ambig_ref.push_back(new_ambig_base);

				// 4. fill emission tables
				double temp;
				for (char j : { 'A', 'C', 'G', 'T' })
				{
					temp = 0;
					for (char k : { 'A', 'C', 'G', 'T' })
					{
						temp += cur_allel_table[k] * (j == k ? 1 - error_rate : error_rate / 3);
					}
					m_E[i][j] = static_cast<T>(log_base((ambig_bases_unequal_weight ? temp : (cur_allel_table[j] ? 1.0 : temp))));
					m_table_of_included_bases[i][j] = (cur_allel_table[j] > 0.0);
				}

				m_E[i]['N'] = static_cast<T>(log_base(1.0));
				m_table_of_included_bases[i]['N'] = true;
				m_allel_freq[i] = cur_allel_table;
			}
			else
			{
				// have no coverage, need to take sacrifical steps to construct
				if (m_emission_tables_initialized == true)
				{
					if ((m_ambig_ref.length() == m_L) && (m_majority_ref.length() == m_L))
					{
						// use old emissions from some other source
						new_majority_ref.push_back(std::tolower(m_majority_ref[i]));
						new_ambig_ref.push_back(std::tolower(m_ambig_ref[i]));
					}
					else
					{
						// error out
						std::cerr << "ERROR: Lengths of old and new emission tables do not match up." << std::endl;
						exit(EXIT_FAILURE);
					}
				}
				else
				{
					// error out
					std::cerr << "ERROR: Emission table not yet initialised and no coverage at position " << i << " of genome." << std::endl;
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	m_majority_ref = std::move(new_majority_ref);
	m_ambig_ref = std::move(new_ambig_ref);
	m_emission_tables_initialized = true;
}

template <typename T>
void reference_genome<T>::init_parameters(
	const background_rates& error_rates,
	const bool ambig_bases_unequal_weight)
{
	T(*cast)
	(const double&) = type_caster<double, T>;

	/* length of HMM */
	m_L = m_allel_freq.size();

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
	if (m_vec_M_to_D_p.size() != m_L - 1)
	{
		std::cerr << "m_vec_M_to_D_p (L = " << m_vec_M_to_D_p.size() << ") needs to have same length as emission matrix (L = " << m_L << ")!" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (m_vec_D_to_D_p.size() != m_L - 1)
	{
		std::cerr << "m_vec_D_to_D_p (L = " << m_vec_D_to_D_p.size() << ") needs to have same length as emission matrix (L = " << m_L << ")!" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<trans_matrix<double>> float_trans_matrix(m_L);

	////////////////////////
	// Position 0 in pHMM //
	////////////////////////
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

	////////////////////////
	// Position 1 in pHMM //
	////////////////////////
	double cur_M_to_D = (m_vec_M_to_D_p[0] ? m_vec_M_to_D_p[0] : error_rates.gap_open);
	double sum = error_rates.insert_open + cur_M_to_D + error_rates.right_clip_open + error_rates.end_prob;
	double M_to_M = 1.0 - sum;
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

	////////////////////////
	// Position i in pHMM //
	////////////////////////
	double cur_D_to_D;
	for (typename std::vector<trans_matrix<double>>::size_type i = 2; i < m_L - 1; ++i)
	{
		cur_M_to_D = (m_vec_M_to_D_p[i - 1] ? m_vec_M_to_D_p[i - 1] : error_rates.gap_open);

		// -> D transitions
		cur_D_to_D = (m_vec_D_to_D_p[i - 1] ? m_vec_D_to_D_p[i - 1] : error_rates.gap_extend);
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

	//////////////////////////
	// Position L-1 in pHMM //
	//////////////////////////
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
	for (typename std::vector<trans_matrix<double>>::size_type i = 0; i < m_L; ++i)
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

	////////////////////////////
	// emission probabilities //
	////////////////////////////
	init_emission_table(
		m_allel_freq,
		error_rates.low_frequency_cutoff,
		error_rates.error_rate,
		ambig_bases_unequal_weight,
		false);

	std::cout << '\t' << "Size of genome:          " << m_L << " nt" << std::endl
			  << std::endl;
	create_index();
}

template <typename T>
void reference_genome<T>::set_parameters(
	const std::string& input_msa,
	const background_rates& error_rates,
	const uint32_t read_lengths,
	const bool ambig_bases_unequal_weight,
	msa_input)
{
	read_length_profile = read_lengths;
	set_parameters(fasta_read<reference_haplotype>(input_msa), error_rates, ambig_bases_unequal_weight);
}

template <typename T>
void reference_genome<T>::set_parameters(
	const std::string& input_hmm_archive,
	const background_rates& error_rates,
	const uint32_t read_lengths,
	const bool ambig_bases_unequal_weight,
	serialized_input)
{
	read_length_profile = read_lengths;
	std::ifstream ifs(input_hmm_archive);
	load_from_file(ifs);

	init_parameters(error_rates, ambig_bases_unequal_weight);
}

template <typename T>
void reference_genome<T>::set_parameters(
	const std::vector<reference_haplotype>& refs,
	const background_rates& error_rates,
	const bool ambig_bases_unequal_weight)
{
	const int32_t L = refs[0].sequence.length();
	const int32_t num_haps = refs.size();

	std::vector<dna_array<double, 5>> E_p(L, { 0.0, 0.0, 0.0, 0.0, 0.0 });
	std::vector<double> M_D_p(L - 1, 0);
	std::vector<double> D_D_p(L - 1, 0);

	// 1.) check that all haplotypes have the same length
	for (int32_t i = 1; i < num_haps; ++i)
	{
		if (L != static_cast<int32_t>(refs[i].sequence.length()))
		{
			std::cerr << "ERROR: Haplotype '" << refs[i].name << "' does not have length L = " << L << "! Aborting." << std::endl;
			exit(EXIT_FAILURE);
		}

		if (refs[i].count <= 0)
		{
			std::cerr << "ERROR: Haplotype '" << refs[i].name << "' has non-positive count " << refs[i].count << "! Aborting." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// 2.) check that all haplotypes have proper DNA letters
	char cur_base;
	for (int32_t i = 0; i < num_haps; ++i)
	{
		for (int32_t j = 0; j < L; ++j)
		{
			cur_base = refs[i].sequence[j];
			if ((wobble_to_ambig_bases.find(cur_base) == wobble_to_ambig_bases.end()) && (cur_base != '-'))
			{
				std::cerr << "ERROR: Unknown base '" << cur_base << "' at position '" << j << "' of '" << refs[i].name << "'." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	// 3.) start producing emission and transition tables
	double sum;
	double MD, sumM;
	double DD, sumD;
	bool onlyGap;

	for (int32_t j = 0; j < L; ++j)
	{
		sum = 0;

		MD = 0;
		sumM = 0;
		DD = 0;
		sumD = 0;

		onlyGap = true;

		// first loop
		for (int32_t i = 0; i < num_haps; ++i)
		{
			if ((refs[i].start <= j) && (j < refs[i].end))
			{
				cur_base = refs[i].sequence[j];
				onlyGap &= (cur_base == '-');

				// allele frequencies
				if (cur_base != '-')
				{
					sum += refs[i].count;

					const std::string& bases = wobble_to_ambig_bases.find(cur_base)->second;
					const uint8_t num_chars = bases.length();
					for (auto k : bases)
					{
						E_p[j][k] += refs[i].count / num_chars;
					}
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
			std::cerr << "ERROR: Position '" << j << "' only has gaps." << std::endl;
			exit(EXIT_FAILURE);
		}

		// renormalize counts
		for (char k : { 'A', 'C', 'G', 'T' })
		{
			E_p[j][k] /= sum;
		}

		if (j < L - 1)
		{
			M_D_p[j] = MD / (sumM ? sumM : 1);
			D_D_p[j] = DD / (sumD ? sumD : 1);
		}
	}

#ifndef NDEBUG
	std::ofstream debug_output("DEBUG_emissions.log");
	for (std::string::size_type i = 0; i < L; ++i)
	{
		debug_output << "Pos " << i << ":\n"
					 << E_p[i] << '\n';
	}
	debug_output.close();
#endif

	set_parameters(std::move(E_p), std::move(M_D_p), std::move(D_D_p), error_rates, ambig_bases_unequal_weight);
}

template <typename T>
void reference_genome<T>::load_from_file(std::istream& ifs)
{
	constexpr auto max_seek = std::numeric_limits<std::streamsize>::max();

	// load length
	ifs.ignore(max_seek, ' ');
	decltype(m_L) L;
	ifs >> L;
	ifs.ignore(max_seek, '\n');

	// ignore "Emission Table" line
	ifs.ignore(max_seek, '\n');

	// ignore "A C G T N" line
	ifs.ignore(max_seek, '\n');

	m_allel_freq.resize(L);
	for (std::size_t i = 0; i < L; ++i)
	{
		ifs >> m_allel_freq[i];
	}
	ifs.ignore(max_seek, '\n');

	// ignore separating line "M -> D Table"
	ifs.ignore(max_seek, '\n');
	m_vec_M_to_D_p.resize(L - 1);
	for (std::size_t i = 0; i < L - 1; ++i)
	{
		ifs >> m_vec_M_to_D_p[i];
	}
	ifs.ignore(max_seek, '\n');

	// ignore separating line "D -> D Table"
	ifs.ignore(max_seek, '\n');
	m_vec_D_to_D_p.resize(L - 1);
	for (std::size_t i = 0; i < L - 1; ++i)
	{
		ifs >> m_vec_D_to_D_p[i];
	}

	// perform final sanity check on stream
	if (ifs.fail() == true)
	{
		std::cerr << "\tERROR during reading of parameters, aborting." << std::endl;
		exit(EXIT_FAILURE);
	}
}

template <typename T>
void reference_genome<T>::save_to_file(std::ostream& ofs)
{
	ofs << std::fixed << std::setprecision(5);
	ofs << "L: " << m_allel_freq.size() << '\n';
	ofs << "Emission Table:\n";
	ofs
		<< std::left << std::setw(5 + 3) << 'A'
		<< std::left << std::setw(5 + 3) << 'C'
		<< std::left << std::setw(5 + 3) << 'G'
		<< std::left << std::setw(5 + 3) << 'T'
		<< std::left << std::setw(5 + 3) << 'N'
		<< '\n';

	for (const auto& i : m_allel_freq)
	{
		ofs << i << '\n';
	}

	ofs << "M -> D Table\n";
	for (const auto& i : m_vec_M_to_D_p)
	{
		ofs << i << '\n';
	}

	ofs << "D -> D Table\n";
	for (const auto& i : m_vec_D_to_D_p)
	{
		ofs << i << '\n';
	}
}

// Misc
template <typename V, typename std::enable_if<std::is_floating_point<V>::value, int>::type = 0>
void inline print_sum(std::ostream&, V v, const bool fail_on_non_summation, const char* suffix = "", const char* prefix = "")
{
	std::cout << prefix << std::right << std::setw(30 - strlen(prefix)) << v << suffix;
	V diff = std::fabs(v - 1.0);
	if ((diff > 1E-4) && (fail_on_non_summation))
	{
		std::cerr << std::setprecision(12) << "ERROR: Failed sum: " << diff << std::endl;
		exit(EXIT_FAILURE);
	}
}

template <typename V, typename std::enable_if<std::is_integral<V>::value, int>::type = 0>
void inline print_sum(std::ostream&, V, const bool, const char* = "", const char* = "")
{
}

template <typename T>
bool reference_genome<T>::display_parameters(std::ostream& output, const bool fail_on_non_summation) const
{
	std::setlocale(LC_ALL, "en_US.UTF-8");
	output
		<< std::setprecision(3) << std::fixed;
	output
		<< "HMM Transition tables for <" << typeid(m_trans_matrix[0].into_match.from_left_clip).name() << ">:" << std::endl;

	constexpr int tr_col_width = 8;
	constexpr int value_col_width = 9;

	/* BEGIN */
	output
		<< std::left << std::setw(tr_col_width) << "BEG"
		<< u8"\u2500\u252C\u2500> "
		<< std::left << std::setw(5 + 2) << "S_L" << ':'
		<< std::right << std::setw(value_col_width) << m_into_left_clip.from_begin << std::endl;
	output
		<< std::left << std::string(tr_col_width, ' ')
		<< u8" \u2514\u2500> "
		<< std::left << std::setw(5 + 2) << "M_j" << ':'
		<< std::right << std::setw(value_col_width) << m_trans_matrix[0].into_match.from_begin << std::endl;
	print_sum<T>(output, m_into_left_clip.from_begin + m_L * m_trans_matrix[0].into_match.from_begin, fail_on_non_summation, "\n", "Sum:");

	/* LEFT CLIP */
	output
		<< std::endl
		<< std::endl
		<< std::left << std::setw(tr_col_width) << "S_L"
		<< u8"\u2500\u252C\u2500> "
		<< std::left << std::setw(5 + 2) << "S_L" << ':'
		<< std::right << std::setw(value_col_width) << m_into_left_clip.from_left_clip << std::endl;
	output
		<< std::left << std::string(tr_col_width, ' ')
		<< u8" \u2514\u2500> "
		<< std::left << std::setw(5 + 2) << "M_j" << ':'
		<< std::right << std::setw(value_col_width) << m_trans_matrix[0].into_match.from_left_clip << std::endl;
	print_sum<T>(output, m_into_left_clip.from_left_clip + m_L * m_trans_matrix[0].into_match.from_left_clip, fail_on_non_summation, "\n", "Sum:");

	/* MAIN STATES */
	for (typename std::vector<trans_matrix<T>>::size_type i = 0; i < m_L - 1; ++i)
	{
		// first row
		output
			<< std::endl
			<< std::endl
			<< std::left << "M_" << std::setw(tr_col_width - 2) << i
			<< u8"\u2500\u252C\u2500> "
			<< std::left << "M_" << std::setw(5) << i + 1 << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_match.from_match;
		output
			<< std::string(2 * tr_col_width, ' ')
			<< std::left << "I_" << std::setw(tr_col_width - 2) << i
			<< u8"\u2500\u252C\u2500> "
			<< std::left << "M_" << std::setw(5) << i + 1 << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_match.from_insertion;

		if (i != 0)
		{
			output
				<< std::string(2 * tr_col_width, ' ')
				<< std::left << "D_" << std::setw(tr_col_width - 2) << i
				<< u8"\u2500\u252C\u2500> "
				<< std::left << "M_" << std::setw(5) << i + 1 << ':'
				<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_match.from_deletion;
		}

		// second row
		output
			<< std::endl
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u251C\u2500> "
			<< std::left << "D_" << std::setw(5) << i + 1 << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_deletion.from_match;
		output
			<< std::string(2 * tr_col_width, ' ')
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u2514\u2500> "
			<< std::left << "I_" << std::setw(5) << i << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i].into_insertion.from_insertion;

		if (i != 0)
		{
			output
				<< std::string(2 * tr_col_width, ' ')
				<< std::left << std::string(tr_col_width, ' ')
				<< u8" \u2514\u2500> "
				<< std::left << "D_" << std::setw(5) << i + 1 << ':'
				<< std::right << std::setw(value_col_width) << m_trans_matrix[i + 1].into_deletion.from_deletion;
		}

		// third row
		output
			<< std::endl
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u251C\u2500> "
			<< std::left << "I_" << std::setw(5) << i << ':'
			<< std::right << std::setw(value_col_width) << m_trans_matrix[i].into_insertion.from_match << std::endl;

		// fourth row
		output
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u251C\u2500> "
			<< std::left << std::setw(5 + 2) << "S_R" << ':'
			<< std::right << std::setw(value_col_width) << m_into_right_clip.from_match << std::endl;

		// fifth row
		output
			<< std::left << std::string(tr_col_width, ' ')
			<< u8" \u2514\u2500> "
			<< std::left << std::setw(5 + 2) << "END" << ':'
			<< std::right << std::setw(value_col_width) << m_into_end.from_match << std::endl;

		print_sum<T>(output,
			m_trans_matrix[i + 1].into_match.from_match + m_trans_matrix[i + 1].into_deletion.from_match + m_trans_matrix[i].into_insertion.from_match + m_into_right_clip.from_match + m_into_end.from_match,
			fail_on_non_summation, "", "Sum:");

		output << std::string(2 * tr_col_width, ' ');
		print_sum<T>(output,
			m_trans_matrix[i + 1].into_match.from_insertion + m_trans_matrix[i].into_insertion.from_insertion,
			fail_on_non_summation);
		if (i != 0)
		{
			output << std::string(2 * tr_col_width, ' ');
			print_sum<T>(output,
				m_trans_matrix[i + 1].into_match.from_deletion + m_trans_matrix[i + 1].into_deletion.from_deletion,
				fail_on_non_summation, "\n");
		}
	}

	/* LAST MATCH STATE */
	output
		<< std::endl
		<< std::endl
		<< std::left << "M_" << std::setw(tr_col_width - 2) << m_L - 1
		<< u8"\u2500\u252C\u2500> "
		<< std::left << std::setw(5 + 2) << "S_R" << ':'
		<< std::right << std::setw(value_col_width) << m_into_right_clip.from_match << std::endl;
	output
		<< std::left << std::string(tr_col_width, ' ')
		<< u8" \u2514\u2500> "
		<< std::left << std::setw(5 + 2) << "END" << ':'
		<< std::right << std::setw(value_col_width) << m_into_end.from_last_match << std::endl;
	print_sum<T>(output, m_into_right_clip.from_match + m_into_end.from_last_match, fail_on_non_summation, "\n", "Sum:");

	/* RIGHT CLIP */
	output
		<< std::endl
		<< std::endl
		<< std::left << std::setw(tr_col_width) << "S_R"
		<< u8"\u2500\u252C\u2500> "
		<< std::left << std::setw(5 + 2) << "S_R" << ':'
		<< std::right << std::setw(value_col_width) << m_into_right_clip.from_right_clip << std::endl;
	output
		<< std::left << std::string(tr_col_width, ' ')
		<< u8" \u2514\u2500> "
		<< std::left << std::setw(5 + 2) << "END" << ':'
		<< std::right << std::setw(value_col_width) << m_into_end.from_right_clip << std::endl;
	print_sum<T>(output, m_into_right_clip.from_right_clip + m_into_end.from_right_clip, fail_on_non_summation, "\n", "Sum:");

	output << std::boolalpha << std::fixed << std::setprecision(5) << std::endl
		   << "Emission tables (log):" << std::endl;
	for (int32_t j = 0; j < m_L; ++j)
	{
		output << std::left << std::setw(7) << j;
		for (char c : { 'A', 'C', 'G', 'T', 'N' })
		{
			output << std::right << std::setw(8) << m_E[j][c] << " (" << m_allel_freq[j][c] << ", " << std::setw(5) << m_table_of_included_bases[j][c] << ")";
		}
		output << std::endl;
	}
	output << std::endl
		   << "Emission tables (uniform):" << std::endl
		   << '\t' << m_uniform_base_e << std::endl;

	return true;
}
}

#endif /* NGSHMMALIGN_REFERENCE_IMPL_HPP */
