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

#pragma once

#include <Rank.hpp>
#include <SimpleSelectHalf.hpp>
#include <SimpleSelectZeroHalf.hpp>
#include <cstdint>
#include <vector>

namespace sux::bits {

using namespace std;
using namespace sux;

/** An implementation of selection and ranking based on the Elias-Fano representation
 * of monotone sequences.
 *
 * Instances of this class can be built using a bit vector or an explicit list of
 * positions for the ones in a vector. In every case, the bit vector or the list
 * are not necessary after construction.
 *
 * @tparam AT a type of memory allocation out of sux::util::AllocType.
 */

template <util::AllocType AT = util::AllocType::MALLOC> class EliasFano : public Rank, public Select {
  private:
	util::Vector<uint64_t, AT> lower_bits, upper_bits;
	SimpleSelectHalf<AT> select_upper;
	SimpleSelectZeroHalf<AT> selectz_upper;
	uint64_t num_bits, num_ones;
	int l;
	int block_size;
	int block_length;
	uint64_t block_size_mask;
	uint64_t lower_l_bits_mask;
	uint64_t ones_step_l;
	uint64_t msbs_step_l;
	uint64_t compressor;
	uint64_t len;

	__inline static void set(util::Vector<uint64_t, AT> &bits, const uint64_t pos) { bits[pos / 64] |= 1ULL << pos % 64; }

	__inline static uint64_t get_bits(util::Vector<uint64_t, AT> &bits, const uint64_t start, const int width) {
		const int start_word = start / 64;
		const int start_bit = start % 64;
		const int total_offset = start_bit + width;
		const uint64_t result = bits[start_word] >> start_bit;
		return (total_offset <= 64 ? result : result | bits[start_word + 1] << (64 - start_bit)) & ((1ULL << width) - 1);
	}

	__inline static void set_bits(util::Vector<uint64_t, AT> &bits, const uint64_t start, const int width, const uint64_t value) {
		const uint64_t start_word = start / 64;
		const uint64_t end_word = (start + width - 1) / 64;
		const uint64_t start_bit = start % 64;

		if (start_word == end_word) {
			bits[start_word] &= ~(((1ULL << width) - 1) << start_bit);
			bits[start_word] |= value << start_bit;
		} else {
			// Here start_bit > 0.
			bits[start_word] &= (1ULL << start_bit) - 1;
			bits[start_word] |= value << start_bit;
			bits[end_word] &= -(1ULL << (width - 64 + start_bit));
			bits[end_word] |= value >> (64 - start_bit);
		}
	}

  public:
	/** Creates a new instance using an
	 *  explicit list of positions for the ones in a bit vector.
	 *
	 *  Note that the list is read only at construction time.
	 *
	 *  In practice this constructor builds an Elias-Fano
	 *  representation of the given list. select(const uint64_t rank) will retrieve
	 *  an element of the list, and rank(const size_t pos) will return how many
	 *  element of the list are smaller than the argument.
	 *
	 * @param ones a list of positions of the ones in a bit vector.
	 * @param num_bits the length (in bits) of the bit vector.
	 */
	EliasFano(const std::vector<uint64_t> ones, const uint64_t num_bits)
	{
		num_ones = ones.size();
		this->num_bits = num_bits;
		l = num_ones == 0 ? 0 : max(0, lambda_safe(num_bits / num_ones));
		this->len = 15;

#ifdef DEBUG
		printf("Number of ones: %lld l: %d\n", num_ones, l);
		printf("Upper bits: %lld\n", num_ones + (num_bits >> l) + 1);
		printf("Lower bits: %lld\n", num_ones * l);
#endif

		const uint64_t lower_bits_mask = (1ULL << l) - 1;
		std::cout << "Lower bit mask: " << lower_bits_mask << std::endl;

		lower_bits.size((num_ones * l + 63) / 64 + 2 * (l == 0));
		upper_bits.size(((num_ones + (num_bits >> l) + 1) + 63) / 64);

		std::cout << "lower_bits size: " << lower_bits.size() << std::endl;
		std::cout << "upper_bits size: " << lower_bits.size() << std::endl;

		for (uint64_t i = 0; i < num_ones; i++) {
			if (l != 0) set_bits(lower_bits, i * l, l, ones[i] & lower_bits_mask);
			set(upper_bits, (ones[i] >> l) + i);
		}

#ifdef DEBUG
		printf("First lower: %016llx %016llx %016llx %016llx\n", lower_bits[0], lower_bits[1], lower_bits[2], lower_bits[3]);
		printf("First upper: %016llx %016llx %016llx %016llx\n", upper_bits[0], upper_bits[1], upper_bits[2], upper_bits[3]);
#endif

		select_upper = SimpleSelectHalf<>(&upper_bits, num_ones + (num_bits >> l) + 1);
		selectz_upper = SimpleSelectZeroHalf<>(&upper_bits, num_ones + (num_bits >> l) + 1);

		block_size = 0;
		do ++block_size;
		while (block_size * l + block_size <= 64 && block_size <= l);
		block_size--;

#ifdef DEBUG
		printf("Block size: %d\n", block_size);
#endif
		block_size_mask = (1ULL << block_size) - 1;
		block_length = block_size * l;

#ifdef PARSEARCH
		ones_step_l = 0;
		for (int i = 0; i < block_size; i++) ones_step_l |= 1ULL << i * l;
		msbs_step_l = ones_step_l << (l - 1);

		compressor = 0;
		for (int i = 0; i < block_size; i++) compressor |= 1ULL << ((l - 1) * i + block_size);
#endif

		lower_l_bits_mask = (1ULL << l) - 1;
		/*
		uint64_t rank_res = rank_par(12);
		std::cout << "----->RANK(19): " << rank_res << std::endl; 
		uint64_t sel_res = select(rank_res);
		std::cout << "----->SELECT(RANK(19)): " << sel_res << std::endl; 
		uint64_t upper = rank_par(sel_res+1);
		std::cout << "----->[" << rank_res << "-" << upper << "]" << std::endl; 
		*/
	}

	uint64_t rank(const size_t k)
	{
		if (num_ones == 0) return 0;
		if (k >= num_bits) return num_ones;
		const uint64_t k_shiftr_l = k >> l;

		int64_t pos = selectz_upper.selectZero(k_shiftr_l);
		uint64_t rank = pos - (k_shiftr_l);

		uint64_t rank_times_l = rank * l;
		const uint64_t k_lower_bits = k & lower_l_bits_mask;

		do {
			rank--;
			rank_times_l -= l;
			pos--;
		} while (pos >= 0 && (upper_bits[pos / 64] & 1ULL << pos % 64) && get_bits(lower_bits, rank_times_l, l) >= k_lower_bits);

		return ++rank;
	}

	uint64_t rank_par(const size_t k) {
		if (num_ones == 0) return 0;
		if (k >= num_bits) return num_ones;
		const uint64_t k_shiftr_l = k >> l;

		const uint64_t k_lower_bits = k & lower_l_bits_mask;

		const uint64_t k_lower_bits_step_l = k_lower_bits * ones_step_l;

		uint64_t pos = selectz_upper.selectZero(k_shiftr_l);
		uint64_t rank = pos - (k_shiftr_l);
		uint64_t rank_times_l = rank * l;

		uint64_t block_upper_bits, block_lower_bits;

		while (rank > block_size) {
			rank -= block_size;
			rank_times_l -= block_length;
			pos -= block_size;
			block_upper_bits = get_bits(upper_bits, pos, block_size);
			block_lower_bits = get_bits(lower_bits, rank_times_l, block_length);

			// printf( "block upper bits: %llx block lower bits: %llx\n", block_upper_bits, block_lower_bits
			// );

			const uint64_t cmp =
				((((block_lower_bits | msbs_step_l) - (k_lower_bits_step_l & ~msbs_step_l)) | (k_lower_bits_step_l ^ block_lower_bits)) ^ (k_lower_bits_step_l & ~block_lower_bits)) & msbs_step_l;

			// printf( "Compare: %016llx compressed: %016llx shifted: %016llx\n", cmp, cmp * compressor, cmp
			// * compressor >> block_size * l );

			const uint64_t cmp_compr = ~(cmp * compressor >> block_length & block_upper_bits) & block_size_mask;

			// printf( "Combined compare: %llx\n", ~t );

			if (cmp_compr) return rank + 1 + lambda_safe(cmp_compr);
		}

		block_upper_bits = get_bits(upper_bits, pos - rank, rank);
		block_lower_bits = get_bits(lower_bits, 0, rank_times_l);

		// printf( "\nTail (%lld bits)...\n", rank );

		// printf( "block upper bits: %llx block lower bits: %llx\n", block_upper_bits, block_lower_bits
		// );

		const uint64_t cmp =
			((((block_lower_bits | msbs_step_l) - (k_lower_bits_step_l & ~msbs_step_l)) | (k_lower_bits_step_l ^ block_lower_bits)) ^ (k_lower_bits_step_l & ~block_lower_bits)) & msbs_step_l;

		// printf( "Compressor: %llx\n", compressor );
		// printf( "Compare: %016llx compressed: %016llx shifted: %016llx\n", cmp, cmp * compressor, cmp *
		// compressor >> block_size * l );

		const uint64_t cmp_compr = ~(cmp * compressor >> block_length & block_upper_bits) & (1ULL << rank) - 1;

		// printf( "Combined compare: %llx\n", ~t );

		return 1 + lambda_safe(cmp_compr);
	}

	/*
	uint64_t rank(const size_t k) {
		if (num_ones == 0) return 0;
		if (k >= num_bits) return num_ones;
#ifdef DEBUG
		printf("Ranking %lld...\n", k);
#endif
		const uint64_t k_shiftr_l = k >> l;
		std::cout << "k: " << k << std::endl;
		std::cout << "l: " << l << std::endl;
		std::cout << "k_shiftr_l: " << k_shiftr_l << std::endl;

#ifndef PARSEARCH
		std::cout << "PARSEARCH" << std::endl;
		int64_t pos = selectz_upper.selectZero(k_shiftr_l);
		int64_t pos2 = selectz_upper.selectZero(k_shiftr_l+1);
		int64_t pos3 = selectz_upper.selectZero(k_shiftr_l+2);
		int64_t pos0 = selectz_upper.selectZero(k_shiftr_l-1);
		std::cout << "--pos0: " << pos0 << std::endl;
		std::cout << "--pos1: " << pos << std::endl;
		std::cout << "--pos2: " << pos2 << std::endl;
		std::cout << "--pos3: " << pos3 << std::endl;
		uint64_t rank = pos - (k_shiftr_l);

		std::cout << "rank: " << rank << std::endl;

#ifdef DEBUG
		printf("Position: %lld rank: %lld\n", pos, rank);
#endif
		uint64_t rank_times_l = rank * l;
		const uint64_t k_lower_bits = k & lower_l_bits_mask;

		do {
			rank--;
			rank_times_l -= l;
			pos--;
		} while (pos >= 0 && (upper_bits[pos / 64] & 1ULL << pos % 64) && get_bits(lower_bits, rank_times_l, l) >= k_lower_bits);

		return ++rank;
#else
		std::cout << "NOT PARSEARCH" << std::endl;

		const uint64_t k_lower_bits = k & lower_l_bits_mask;

#ifdef DEBUG
		printf("k: %llx lower %d : %llx\n", k, l, k_lower_bits);
#endif

		const uint64_t k_lower_bits_step_l = k_lower_bits * ones_step_l;

		uint64_t pos = selectz_upper.selectZero(k_shiftr_l);
		uint64_t rank = pos - (k_shiftr_l);
		uint64_t rank_times_l = rank * l;

#ifdef DEBUG
		printf("pos: %lld rank: %lld\n", pos, rank);
#endif

		uint64_t block_upper_bits, block_lower_bits;

		while (rank > block_size) {
			rank -= block_size;
			rank_times_l -= block_length;
			pos -= block_size;
			block_upper_bits = get_bits(upper_bits, pos, block_size);
			block_lower_bits = get_bits(lower_bits, rank_times_l, block_length);

			// printf( "block upper bits: %llx block lower bits: %llx\n", block_upper_bits, block_lower_bits
			// );

			const uint64_t cmp =
				((((block_lower_bits | msbs_step_l) - (k_lower_bits_step_l & ~msbs_step_l)) | (k_lower_bits_step_l ^ block_lower_bits)) ^ (k_lower_bits_step_l & ~block_lower_bits)) & msbs_step_l;

			// printf( "Compare: %016llx compressed: %016llx shifted: %016llx\n", cmp, cmp * compressor, cmp
			// * compressor >> block_size * l );

			const uint64_t cmp_compr = ~(cmp * compressor >> block_length & block_upper_bits) & block_size_mask;

			// printf( "Combined compare: %llx\n", ~t );

			if (cmp_compr) return rank + 1 + lambda_safe(cmp_compr);
		}

		block_upper_bits = get_bits(upper_bits, pos - rank, rank);
		block_lower_bits = get_bits(lower_bits, 0, rank_times_l);

		// printf( "\nTail (%lld bits)...\n", rank );

		// printf( "block upper bits: %llx block lower bits: %llx\n", block_upper_bits, block_lower_bits
		// );

		const uint64_t cmp =
			((((block_lower_bits | msbs_step_l) - (k_lower_bits_step_l & ~msbs_step_l)) | (k_lower_bits_step_l ^ block_lower_bits)) ^ (k_lower_bits_step_l & ~block_lower_bits)) & msbs_step_l;

		// printf( "Compressor: %llx\n", compressor );
		// printf( "Compare: %016llx compressed: %016llx shifted: %016llx\n", cmp, cmp * compressor, cmp *
		// compressor >> block_size * l );

		const uint64_t cmp_compr = ~(cmp * compressor >> block_length & block_upper_bits) & (1ULL << rank) - 1;

		// printf( "Combined compare: %llx\n", ~t );

		return 1 + lambda_safe(cmp_compr);

#endif
	}
	*/

	/*
	std::tuple<uint64_t,uint64_t,bool> 
		binary_search_lower_bound(const std::string& P, uint64_t b, uint64_t e) const
	{
		uint64_t plen = e-b;
		uint64_t to_match = std::min(static_cast<uint_t>(this->len),plen);
		uint64_t r, r_;
		{ // Run Elias-Fano search
			uint64_t z; 
			bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&z), P, len, e-to_match, to_match);
			// search for the bitpacked pattern
			r = rank(z);
			uint64_t s = select(r);

		}
	
		{ // Run the Elias-Fano search
			// bitpack first to_match characters
			uint_t search = 0;
			this->bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
				                                                   pend-to_match, to_match);
			// search for the bitpacked pattern
			       r = ef.rank1(search);    
			uint_t s = ef.select1(r);
			uint_t mlen = (__builtin_clz(search ^ s)-((sizeof(uint_t)*8)-(this->len*2)))/2;

			if( mlen < to_match )
				return std::make_tuple(-1,0,true);
		
			if(to_match == this->len)
			{ 
				r_ = r+1;
			}
			else
			{
				uint_t search_ = search | (1 << (this->len-to_match)*2)-1;
				this->bitpack_uint_DNA(reinterpret_cast<uint8_t*>(&search_), pattern, this->len, 
					                                                   pend-to_match, to_match);
				r_ = ef.rank1(search_)-1;
			}
			low  = bv.select(r); high = bv.select(r_);

			if(this->get_sample(low)+1 < this->len)
			{
				r++;
				s = ef.select1(r);
				mlen = (__builtin_clz(search ^ s)-((sizeof(uint_t)*8)-(this->len*2)))/2;

				if( mlen < to_match )
					return std::make_tuple(-1,0,0);

				low  = bv.select(r);
				if(r == r_){ r_++; high = bv.select(r_); }
			}
		}
		{ // Run the binary search
			// stop if first pattern character doesn't occur in the text
			if((high - low) > 0)
			{ 
				high--;
				lcp = O->LCS(pattern,pend-1,this->get_sample(high)); 
				mid = (low+high)/2;
			}
			else
				return std::make_tuple(-1,0,true);

			while( low < high )
			{	
				auto j = O->LCS_char(pattern,pend-1,this->get_sample(mid)); 
		
				if((j.first != plen) and (j.second < pattern[pend-j.first-1]))    
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
		}

		return std::make_tuple(this->get_sample(low),lcp,(lcp != plen));
		
		return std::make_tuple(0,0,0);
	}*/

	size_t select(const uint64_t rank)
	{
		return (select_upper.select(rank) - rank) << l | get_bits(lower_bits, rank * l, l);
	}

	uint64_t select(const uint64_t rank, uint64_t *const next) {
		uint64_t s, t;
		s = select_upper.select(rank, &t) - rank;
		t -= rank + 1;

		const uint64_t position = rank * l;
		*next = t << l | get_bits(lower_bits, position + l, l);
		return s << l | get_bits(lower_bits, position, l);
	}

	/** Returns the size in bits of the underlying bit vector. */
	size_t size() const { return num_bits; }

	/** Returns an estimate of the size in bits of this structure. */
	uint64_t bitCount() {
		return upper_bits.bitCount() - sizeof(upper_bits) * 8 + lower_bits.bitCount() - sizeof(lower_bits) * 8 + select_upper.bitCount() - sizeof(select_upper) * 8 + selectz_upper.bitCount() -
			   sizeof(selectz_upper) * 8 + sizeof(*this) * 8;
	}
};

} // namespace sux::bits