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

/* Modifications Copyright (C) 2025 Davide Cenzato */

#pragma once

#include <Rank.hpp>
#include <SimpleSelectHalf.hpp>
#include <SimpleSelectZeroHalf.hpp>
#include <cstdint>
#include <vector>
#include <RLZ_DNA_sux.hpp>

namespace sux::bits {

using namespace std;
using namespace sux;

/** An implementation of selection and ranking based on the Elias-Fano representation
 *  of a monotone sequence of integers interleaved by a fixed amount of bites.
 *
 * @tparam AT a type of memory allocation out of sux::util::AllocType.
 */

template <util::AllocType AT = util::AllocType::MALLOC> class InterleavedEliasFano : public Rank, public Select {
  private:
	util::Vector<uint64_t, AT> lower_bits, upper_bits;
	SimpleSelectHalf<AT> select_upper;
	SimpleSelectZeroHalf<AT> selectz_upper;
	uint64_t u, n;
	int l, w;
	uint64_t lower_l_bits_mask;
	uint8_t u_width;

	__inline static void set(util::Vector<uint64_t, AT> &bits, const uint64_t pos) { bits[pos / 64] |= 1ULL << pos % 64; }

	__inline static uint64_t get_bits(const util::Vector<uint64_t, AT> &bits, const uint64_t start, const int width) {
		const int start_word = start / 64;
		const int start_bit = start % 64;
		const int total_offset = start_bit + width;
		const uint64_t result = bits[start_word] >> start_bit;
		return (total_offset <= 64 ? result : result | bits[start_word + 1] << (64 - start_bit)) & ((1ULL << width) - 1);
	}

	__inline static void 
	set_bits(util::Vector<uint64_t, AT> &bits, const uint64_t start, const int width, const uint64_t value)
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

  public:
  /* empty constructor */
  InterleavedEliasFano(){ }

	/** Creates a Elias-Fano compressed sorted multimap using an 
	 *  explicit list of (key, value) pairs where the keys are given
	 *  as an monotonically increasing list of integers.
	 *
	 *  Note that the list is read only at construction time.
	 *
	 *  In practice this function builds an Elias-Fano compressed sorted
	 *  dictionary admitting duplicate keys.
	 *  In short, select(const uint64_t rank) will retrieve the position ith
	 *  largest key, and rank(const size_t pos) will return how many keys in 
	 *  the dictionary are smaller than the argument.
	 *
	 * @param keys_values a list of (key, value) pairs in monotonically increasing order.
	 * @param universe_size size of the largest key that can be represented.
	 * @param width_values number of bits need to represent each values.
	 */

 	void build(const std::vector<std::pair<uint64_t,uint64_t>>& keys_values,
 		         const uint64_t universe_size, const uint8_t  values_width)
	{
		this->n = keys_values.size();
		this->u = universe_size;
		this->l = n == 0 ? 0 : max(0, lambda_safe(u / n));
		this->w = values_width;

    #ifdef DEBUG
      std::cout << "Universe size: " << u << std::endl;
			std::cout << "Number of integers: " << n << " l: " << l << std::endl;
			std::cout << "Upper bits: " << n + (u >> l) + 1 << std::endl;
			std::cout << "Lower bits: " << n * (l + w) << std::endl;
		#endif

		const uint64_t lower_bits_mask = (1ULL << l) - 1;

		// parchè + 2 * (l == 0) ? 
		lower_bits.size(((n * (l + w)) + 63) / 64 );
		upper_bits.size(((n + (u >> l) + 1) + 63) / 64);

		for (uint64_t i = 0; i < n; ++i)
		{
			if (l != 0) set_bits(lower_bits, i * (l + w), l, keys_values[i].first & lower_bits_mask);
			set_bits(lower_bits,(i * (l + w))+l, w, keys_values[i].second);
			set(upper_bits, (keys_values[i].first >> l) + i);
		}

		select_upper = SimpleSelectHalf<>(&upper_bits, n + (u >> l) + 1);
		selectz_upper = SimpleSelectZeroHalf<>(&upper_bits, n + (u >> l) + 1);

		this->lower_l_bits_mask = (1ULL << l) - 1;
	}

	uint64_t rank1(const size_t k) const
	{
		if (n == 0) return 0;
		if (k >= u) return n;
		#ifdef DEBUG
				printf("Ranking %lld...\n", k);
		#endif
		const uint64_t k_shiftr_l = k >> l;

		int64_t pos = selectz_upper.selectZero(k_shiftr_l);
		uint64_t rank = pos - (k_shiftr_l);

		#ifdef DEBUG
				printf("Position: %lld rank: %lld\n", pos, rank);
		#endif

		uint64_t rank_times_l = rank * (l + w);
		const uint64_t k_lower_bits = k & lower_l_bits_mask;

		do {
			rank--;
			rank_times_l -= (l + w);
			pos--;
		} while (pos >= 0 && (upper_bits[pos / 64] & 1ULL << pos % 64) && 
			                   get_bits(lower_bits, rank_times_l, l) >= k_lower_bits);

		return ++rank;
	}

	size_t select1(const uint64_t rank) const
	{
		#ifdef DEBUG
			std::cout << "Selecting " << rank << "..." << std::endl;
			std::cout << "Returning " << (select_upper.select(rank) - rank) << l | get_bits(lower_bits, rank * l, l) <<
			             " = " << select_upper.select(rank) - rank, l << " | " << get_bits(lower_bits, rank * l, l) << std::endl;
		#endif

		return (select_upper.select(rank) - rank) << l | get_bits(lower_bits, rank * (l + w), l);
	}

	size_t select1_value(const uint64_t rank, uint64_t& value) const
	{
		uint64_t interleaved = get_bits(lower_bits, rank * (l + w), (l + w));

		value = interleaved >> l;
		return (select_upper.select(rank) - rank) << l | ( ((1ULL << l) - 1) & interleaved );
	}

	size_t get_sample(const uint64_t i) const { return get_bits(lower_bits, (i * (l + w))+l, w); }

	// return the successor of i and its interleaved value
	std::pair<uint64_t,uint64_t> successor_value(uint64_t i) const
	{
		std::pair<uint64_t,uint64_t> res;
		res.first = select1_value(rank1(i),res.second);
		
		return res;
	}

	int64_t lower_bound(uint64_t key, uint8_t key_width) const
	{
		uint64_t r, s, val;

		r = rank1(key);
		s = select1_value(r, val) ^ key;

		uint8_t mbits = (s == 0) ? key_width : (__builtin_clzll(s) & ~1) - (64 - u_width);
		auto neg_val = ~val + 1;
		val = (key_width <= mbits) ? val : neg_val;

		return val;
	}

	int64_t lower_bound_offset(uint64_t key, uint8_t key_width, uint64_t offset) const
	{
		uint64_t r, s, val;

		r = rank1(key) + offset;
		s = select1_value(r, val) ^ key;

		uint8_t mbits = (s == 0) ? key_width : (__builtin_clzll(s) & ~1) - (64 - u_width);
		val = (key_width <= mbits) ? val : (~val + 1);

		return val;
	}

	std::tuple<int64_t,uint64_t,uint64_t>
	lower_upper_bound_exact(uint64_t key) const
	{
		uint64_t r, r_, s, val;

		r = rank1(key);    
		s = select1_value(r, val);

		if(s != key) return std::make_tuple(-1,0,0);

		r_ = rank1(s+1);

		if(val+1 < u_width/2)
		{
			r++;
			s = select1_value(r, val);

			if(s != key) return std::make_tuple(-1,0,0);

			if(r == r_){ r_++; }
		}

		return std::make_tuple(r,r_,val);
	}

	template<class RAoracle = RLZ_DNA_sux<>>
	std::pair<int64_t,int64_t>
	binary_search_text_oracle(const std::string& P, uint64_t b, uint64_t e,
		                        uint64_t low, uint64_t high, uint8_t lcs_width,
		                                                      RAoracle* oracle) const
	{
		uint64_t mid = (low + high)/2,
		         lcp = oracle->LCS(P, e-1, get_sample(high) >> lcs_width);

		while( low < high )
		{	
			auto j = oracle->LCS_char(P, e-1, get_sample(mid) >> lcs_width); 
	
			if((j.first != (e - b)) and (j.second < P[e-j.first-1]))    
			{
				low = mid+1;
			}
			else
			{
				high = mid;
				lcp = j.first;
			}
 			
			mid = (low+high)/2;
		}

		return std::make_pair(get_sample(low) >> lcs_width,lcp);
	}
	
	/** Returns the number of integers reprenseted in this structure. */
	size_t size() const { return n; }

	/** Returns the universe size. */
	size_t universe_size() const { return u; }

	/** Returns an estimate of the size in bits of this structure. */
	uint64_t bitCount() {
		return upper_bits.bitCount() - sizeof(upper_bits) * 8 + lower_bits.bitCount() - sizeof(lower_bits) * 8 + select_upper.bitCount() - sizeof(select_upper) * 8 + selectz_upper.bitCount() -
			   sizeof(selectz_upper) * 8 + sizeof(*this) * 8;
	}

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
		w_bytes += select_upper.serialize(out);
		w_bytes += selectz_upper.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in)
	{
		in.read((char*)&u, sizeof(u));
		in.read((char*)&n, sizeof(n));
		in.read((char*)&l, sizeof(l));
		in.read((char*)&w, sizeof(w));
		this->lower_l_bits_mask = (1ULL << this->l) - 1;
		this->u_width = (63 - __builtin_clzll(u));

		lower_bits.load(in);
		upper_bits.load(in);
		select_upper.load(in,&upper_bits);
		selectz_upper.load(in,&upper_bits);
	}
};

} // namespace sux::bits