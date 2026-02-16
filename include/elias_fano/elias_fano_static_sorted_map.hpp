/*
 * Sux: Succinct data structures
 *
 * Copyright (C) 2007-2020 Sebastiano Vigna
 *
 *  This library is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * Under Section 7 of GPL version 3, you are granted additional permissions
 * described in the GCC Runtime Library Exception, version 3.1, as published by
 * the Free Software Foundation.
 *
 * You should have received a copy of the GNU General Public License and a copy of
 * the GCC Runtime Library Exception along with this program; see the files
 * COPYING3 and COPYING.RUNTIME respectively.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

/* Code adapted by Davide Cenzato 2026 */

#pragma once

#include <Rank.hpp>
#include <SimpleSelectHalf.hpp>
#include <SimpleSelectZeroHalf.hpp>
#include <cstdint>
#include <vector>

namespace stpd::ef {

using namespace std;
using namespace sux;

/** An implementation of a static sorted multimap based on Elias–Fano encoding,
 *  where a monotone sequence of integer keys (keys may be repeated; uniqueness
 *  is determined by the combination of the key and its position in the sequences) 
 *  is mapped to integer values.
 *
 * @tparam AT a type of memory allocation out of sux::util::AllocType.
 */
template <util::AllocType AT = util::AllocType::MALLOC>
class EliasFanoStaticSortedMap {

  private:

  /* Elias-Fano data structures */
	util::Vector<uint64_t, AT> lower_bits, upper_bits;
	bits::SimpleSelectHalf<AT> select_upper;
	bits::SimpleSelectZeroHalf<AT> selectz_upper;
	/* Ausiliar variables */
	uint64_t u, n;
	int l, w;
	uint64_t lower_l_bits_mask;
	uint8_t u_width; 

  public:

	/* Empty constructor */ 
	EliasFanoStaticSortedMap(){ } 

	/** Creates a Elias-Fano compressed sorted multimap using an 
	 *  explicit list of (key, value) pairs where the keys are given
	 *  as an monotonically increasing list of integers.
	 *
	 *  Note that the list is read only at construction time.
	 *
	 *  In practice this function builds an Elias-Fano compressed sorted
	 *  dictionary admitting duplicate keys.
	 *  In short, select(const uint64_t i) will retrieve the position of the
	 *  ith key in the list, and rank(const size_t val) will return how many  
	 *  keys in the dictionary are smaller than val.
	 *
	 * @param keys_values a list of (key, value) pairs in monotonically increasing order.
	 * @param universe_size size of the largest key that can be represented.
	 * @param largest_value size of the largest value that can be represented.
	 */
 	void build(const std::vector<std::pair<uint64_t,uint64_t>>& keys_values,
 		         								   								 uint64_t largest_key   = 0,
 		                                           uint64_t largest_value = 0)
	{
		if(not largest_key or not largest_value)
		{
			for (const auto& i : keys_values)
				{ 
					largest_key = std::max(largest_key,i.first);
					largest_value = std::max(largest_value,i.second);
				}
		}
		assert(largest_key > 0);

		this->n = keys_values.size();
		this->u = largest_key+1;
		this->w = lambda_safe(largest_value)+1;
		this->l = n == 0 ? 0 : max(0, lambda_safe(u / n));

		#ifdef DEBUG
		  std::cout << "Universe size (u): " << u << std::endl;
		  std::cout << "Values width (w): " << w << std::endl;
			std::cout << "Number of integers (n): " << n << std::endl;
			std::cout << "log(u/n): " << l << std::endl;
			std::cout << "No. upper bits: " << n + (u >> l) + 1 << std::endl;
			std::cout << "No. lower bits: " << n * (l + w) << std::endl;
		#endif

		const uint64_t lower_bits_mask = (1ULL << l) - 1;

		lower_bits.size(((n * (l + w)) + 63) / 64 );
		upper_bits.size(((n + (u >> l) + 1) + 63) / 64);

		for (uint64_t i = 0; i < n; ++i) 
		{
			if (l != 0) set_lower(lower_bits, i * (l + w), l, keys_values[i].first & lower_bits_mask);
			set_lower(lower_bits,(i * (l + w))+l, w, keys_values[i].second);
			set_higher(upper_bits, (keys_values[i].first >> l) + i);
		}

		select_upper = bits::SimpleSelectHalf<>(&upper_bits, n + (u >> l) + 1);
		selectz_upper = bits::SimpleSelectZeroHalf<>(&upper_bits, n + (u >> l) + 1);

		this->lower_l_bits_mask = (1ULL << l) - 1;
		this->u_width = (63 - __builtin_clzll(u));
	}

	/* return the i-th key in key-sorted order (duplicates allowed) */
	uint64_t select_key(const uint64_t rank) const
	{
		#ifdef DEBUG
			std::cout << "Selecting " << rank << "..." << std::endl;
			std::cout << "Returning " << 
						 (select_upper.select(rank) - rank) << l | get_bits(lower_bits, rank * l, l) << " = " <<
			             (select_upper.select(rank) - rank << l) << " | " << get_bits(lower_bits, rank * l, l) << std::endl;
		#endif 

		return (select_upper.select(rank) - rank) << l | get_bits(lower_bits, rank * (l + w), l);
	}

	/* return the i-th key-value pair in key-sorted order (duplicates allowed) */
	std::pair<uint64_t,uint64_t> select(const uint64_t rank) const
	{ 
		#ifdef DEBUG 
			std::cout << "Selecting " << rank << "..." << std::endl;
			std::cout << "Returning " << 
						 (get_bits(lower_bits, rank * (l + w), (l + w)) >> l) << " = " <<
						  get_bits(lower_bits, rank * (l + w), (l + w)) << " >> " << l << " and " << std::endl <<
						 (select_upper.select(rank) - rank) << l | get_bits(lower_bits, rank * l, l) << " = "  <<
			       (select_upper.select(rank) - rank << l) << " | " << get_bits(lower_bits, rank * l, l) << std::endl;
		#endif 

		//uint64_t key = (select_upper.select(rank) - rank) << l | (((1ULL << l) - 1) & interleaved);
		uint64_t key = (select_upper.select(rank) - rank) << l | get_bits(lower_bits, rank * (l + w), l);
		//uint64_t value = interleaved >> l;
		uint64_t value = get_bits(lower_bits, (rank * (l + w)) + l, w);

		return std::make_pair(key, value);
	}

	/* return the number of keys strictly smaller than the given key */
	uint64_t rank(const size_t val) const
	{
		if (n == 0) return 0;
		if (val >= u) return n;
		#ifdef DEBUG
				printf("Ranking %lld...\n", val);
		#endif
		const uint64_t k_shiftr_l = val >> l;

		int64_t pos = selectz_upper.selectZero(k_shiftr_l);
		uint64_t rank = pos - (k_shiftr_l);

		#ifdef DEBUG
				printf("Position: %lld rank: %lld\n", pos, rank);
		#endif

		uint64_t rank_times_l = rank * (l + w);
		const uint64_t k_lower_bits = val & lower_l_bits_mask;

		do {
			rank--;
			rank_times_l -= (l + w);
			pos--;
		} while (pos >= 0 && (upper_bits[pos / 64] & 1ULL << pos % 64) && 
			                   get_bits(lower_bits, rank_times_l, l) >= k_lower_bits);

		return ++rank;
	}

	/* return the first occurrence of a key ≥ the given key in sorted order */
	std::pair<uint64_t,uint64_t> successor(const size_t key) const
	{
		assert(n > 0);
		assert(key < u);
		const uint64_t k_shiftr_l = key >> l;

		int64_t pos = selectz_upper.selectZero(k_shiftr_l);
		uint64_t rank = pos - (k_shiftr_l);

		uint64_t rank_times_l = rank * (l + w);
		const uint64_t val_lower_bits = key & lower_l_bits_mask;

		do {
			rank--;
			rank_times_l -= (l + w);
			pos--;
		}
		while (pos >= 0 && (upper_bits[pos / 64] & 1ULL << pos % 64) && 
			     get_bits(lower_bits, rank_times_l, l) >= val_lower_bits);

		pos++;
		uint64_t x = upper_bits[pos / 64] >> (pos % 64);
		uint64_t trailing_zeros = 0;
		while(x == 0)
		{
			trailing_zeros += 64 - (pos % 64);

			pos += 64 - (pos % 64);
			x = upper_bits[pos / 64] >> (pos % 64);
		}
		trailing_zeros += __builtin_ctzll(x); 

		/*
		return std::make_pair(
					 (k_shiftr_l + trailing_zeros) << l | ((1ULL << l) - 1) & value,
					 value >> l
					 );*/
		return std::make_pair(
					 (k_shiftr_l + trailing_zeros) << l | get_bits(lower_bits, (rank + 1) * (l + w), l),
					 get_bits(lower_bits, ((rank + 1) * (l + w)) + l, w)
					 );
	}

	std::tuple<uint64_t,uint64_t,uint64_t> 
	key_range(const size_t key) const
	{
		assert(n > 0);
		assert(key < u);
		const uint64_t k_shiftr_l = key >> l;

		int64_t pos = selectz_upper.selectZero(k_shiftr_l);
		uint64_t rank = pos - (k_shiftr_l);
		uint64_t rank_times_l = rank * (l + w);

		const uint64_t val_lower_bits = key & lower_l_bits_mask;

		do {
			rank--;
			rank_times_l -= (l + w);
			pos--;
		}
		while (pos >= 0 && (upper_bits[pos / 64] & 1ULL << pos % 64) && 
			     get_bits(lower_bits, rank_times_l, l) > val_lower_bits);
		uint64_t upper_bound = rank + 1;

		while (pos >= 0 && (upper_bits[pos / 64] & 1ULL << pos % 64) && 
			     get_bits(lower_bits, rank_times_l, l) == val_lower_bits)
		{
			rank--;
			rank_times_l -= (l + w);
			pos--;
		}
		uint64_t lower_bound = rank + 1;
	  //uint64_t value = get_bits(lower_bits, (rank + 1) * (l + w), (l + w));

		return std::make_tuple(lower_bound, upper_bound,
					 get_bits(lower_bits, ((rank + 1) * (l + w)) + l, w)
					 //get_bits(lower_bits, (rank + 1) * (l + w), (l + w)) >> l
				   );
	}

	size_t inline get_value(const uint64_t i) const { return get_bits(lower_bits, (i * (l + w))+l, w); }

	uint8_t universe_width() const { return u_width; } 

	size_t serialize(std::ostream& out)
	{
		size_t w_bytes = 0;

		out.write((char*)&u, sizeof(u));
		out.write((char*)&n, sizeof(n));
		out.write((char*)&l, sizeof(l));
		out.write((char*)&w, sizeof(w));

		w_bytes += sizeof(u) + sizeof(n) + sizeof(l) + sizeof(w);

		w_bytes += lower_bits.serialize(out);
		w_bytes += upper_bits.serialize(out);
		w_bytes += selectz_upper.serialize(out);
		w_bytes += select_upper.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, bool successor_only = false)
	{
		in.read((char*)&u, sizeof(u));
		in.read((char*)&n, sizeof(n));
		in.read((char*)&l, sizeof(l));
		in.read((char*)&w, sizeof(w));
		this->lower_l_bits_mask = (1ULL << this->l) - 1;
		this->u_width = (63 - __builtin_clzll(u));

		lower_bits.load(in);
		upper_bits.load(in);
		selectz_upper.load(in,&upper_bits);
		if(not successor_only)
			select_upper.load(in,&upper_bits);
	}

  private:

	inline static void 
	set_higher(util::Vector<uint64_t, AT> &bits, const uint64_t pos) { bits[pos / 64] |= 1ULL << pos % 64; }

	inline static void 
	set_lower(util::Vector<uint64_t, AT> &bits, const uint64_t start, const int width, const uint64_t value)
	{
		const uint64_t start_word = start / 64;
		const uint64_t end_word = (start + width - 1) / 64;
		const uint64_t start_bit = start % 64;

		if (start_word == end_word)
		{
			bits[start_word] &= ~(((1ULL << width) - 1) << start_bit);
			bits[start_word] |= value << start_bit;
		} 
		else
		{
			// Here start_bit > 0.
			bits[start_word] &= (1ULL << start_bit) - 1;
			bits[start_word] |= value << start_bit;
			bits[end_word] &= -(1ULL << (width - 64 + start_bit));
			bits[end_word] |= value >> (64 - start_bit);
		}
	}

	inline static uint64_t
	get_bits(const util::Vector<uint64_t, AT> &bits, const uint64_t start, const int width)
	{
		const uint64_t start_word = start / 64;
		const uint64_t start_bit = start % 64;
		const uint64_t total_offset = start_bit + width;
		const uint64_t result = bits[start_word] >> start_bit;
		return (total_offset <= 64 ? result : result | bits[start_word + 1] << (64 - start_bit)) & ((1ULL << width) - 1);
	}
};

};