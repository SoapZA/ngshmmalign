#ifndef INDEX_IMPL_HPP
#define INDEX_IMPL_HPP

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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace
{

// kmer-based lookup
template <typename T>
void reference_genome<T>::create_index(const uint16_t desired_kmer_length)
{
	//m_kmer_length = desired_kmer_length;
	constexpr int32_t max_kmers = 10000000;

	auto expand_ambig_sequence = [&](const uint32_t POS, const uint32_t length)
	{
		std::stack<char> S;
		char v;
		std::string cur_string;
		uint16_t cur_length = 0;

		for (const auto c : { 'A', 'C', 'G', 'T' })
		{
			if (m_table_of_included_bases[POS][c])
			{
				S.push(c);
			}
		}

		while (!S.empty())
		{
			v = S.top();
			S.pop();

			if (v == 'z')
			{
				// going back in string
				cur_string.pop_back();
				--cur_length;
			}
			else
			{
				cur_string.push_back(v);
				++cur_length;
				S.push('z');

				if (cur_length == length)
				{
					// got a full-length kmer
					auto it = m_kmer_index.find(cur_string);

					if (it == m_kmer_index.end())
					{
						// new kmer not in hash map yet
						m_kmer_index.emplace(cur_string, std::vector<uint32_t>{ POS });
					}
					else
					{
						// kmer already in map
						it->second.push_back(POS);
					}
				}
				else
				{
					// no full length kmer yet
					for (const auto c : { 'A', 'C', 'G', 'T' })
					{
						if (m_table_of_included_bases[POS + cur_length][c])
						{
							S.push(c);
						}
					}
				}
			}
		}
	};

	bool too_large;
	for (m_kmer_length = desired_kmer_length; m_kmer_length >= 10; --m_kmer_length)
	{
		std::cout << "   Building k-mer index, k = " << m_kmer_length;
		m_kmer_index.clear();

		too_large = false;
		for (std::size_t i = 0; i < m_L - m_kmer_length + 1; ++i)
		{
			expand_ambig_sequence(i, m_kmer_length);

			if (m_kmer_index.size() > max_kmers)
			{
				too_large = true;
				break;
			}
		}

		if (too_large == true)
		{
			std::cout << " -> too large\n";
		}
		else
		{
			std::cout << " -> " << m_kmer_index.size() << '\n';
			break;
		}
	}
}

template <typename T>
typename reference_genome<T>::index_stat reference_genome<T>::find_pos(const boost::string_ref& query) const
{
	boost::accumulators::accumulator_set<int64_t, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> acc;
	hash_map_type::const_iterator it;

	const std::size_t max_len = query.length() - m_kmer_length + 1;
	uint32_t num_samples = 0;

	for (std::size_t i = 0; i < max_len; ++i)
	{
		boost::string_ref cur_kmer = query.substr(i, m_kmer_length);
		it = m_kmer_index.find(cur_kmer, m_kmer_index.hash_function(), m_kmer_index.key_eq());

		if (it != m_kmer_index.end())
		{
			for (const auto& j : it->second)
			{
				++num_samples;
				acc(static_cast<int64_t>(j - i));
			}
		}
	}

	int32_t POS;
	int32_t sd;
	if (num_samples)
	{
		// add statistics
		POS = boost::accumulators::mean(acc);
		sd = std::sqrt(boost::accumulators::variance(acc));
	}
	else
	{
		// didn't find anything, likely wrong orientation
		POS = 0;
		sd = 0;
	}

	return { POS, sd, num_samples };
}
}

#endif /* INDEX_IMPL_HPP */