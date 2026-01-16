#ifndef UNI_PART_MOVE_R_HPP_
#define UNI_PART_MOVE_R_HPP_

#include <elias_fano_sux.hpp>
#include <common.hpp>

#include "phi-breakblocks.hpp"

namespace stpd
{
template<class bp_vector = std::vector<uint64_t>, class bitvector = sux::bits::EliasFano<>>
class uni_part_move_r
{
private:
	// array containing the move blocks
	bp_vector blocks;
	// bitvector containing the superblock boundaries
	bitvector borders;

	uint8_t uni_w;
	uint64_t  L, N, U;
	static constexpr int M = 7;
	static constexpr int bit_width = sizeof(typename bp_vector::value_type) * 8;

public:	
	uni_part_move_r(){ } 

	uni_part_move_r(       std::vector<pair_t>&  phi_blocks, 
		        	 const safe_t                universe_size,
		        	 const safe_t                superblock_length = 0,
		        	 const bool_t                verbose = false)
	{ build(phi_blocks, universe_size, superblock_length, verbose); }

	/* expected number of blocks: (r+u/L)*(d+1/d-1) */
	void build(       std::vector<pair_t>&  phi_blocks, 
		        const safe_t                universe_size,
		        const safe_t                superblock_length = 0,
		        const bool_t                verbose = false)
	{
		/* variables initialization */
		U = universe_size;
		safe_t r = phi_blocks.size();
		uint8_t block_w = 64 - __builtin_clzll(r);
		/* phase 1) Break move-r blocks according to the universe size */
		std::vector<pair_t> phi_blocks_u;
		{
			if(superblock_length)
			{
				uni_w  = 64 - __builtin_clzll(superblock_length);
				L = 1ULL << uni_w;
				if(bit_width < (block_w + 2*uni_w + 1))
				{
					std::cerr << "	### Error! The current block vector of width " <<
					             bit_width << " cannot contain the move blocks" <<
					             " correctly. Please select a smaller superblock length or" <<
					             " let this function compute it automatically..." << std::endl;
					             exit(1);
				}
			}
			else
			{
				/* compute the superblock length automatically */
				int_t available_bits = bit_width - 1 - block_w;
				while(true)
				{
					uni_w = available_bits/2;
					L = 1ULL << uni_w;

					safe_t expected_r = (r + U/L) * ((M+1)/double(M-1));
					      expected_r += expected_r/M + M;

					if( (uni_w*2 + 1 + (64 - __builtin_clzll(expected_r))) <= bit_width )
						break; 
					available_bits -= 2;
				}
			}
			/* break the blocks according to the superblock length */
			phi_blocks_u.reserve(r + U/L);
			{
				size_t sum = 0;
				int64_t clen = L;
				for(const auto& i : phi_blocks)
				{
					size_t curr_img_pos = i.first;
					size_t curr_block_len = i.second;
					while(curr_block_len >= clen)
					{
						phi_blocks_u.push_back(std::make_pair(curr_img_pos,clen));
						curr_img_pos += clen;
						curr_block_len -= clen;
						clen = L;
					}
					if(curr_block_len > 0)
					{
						phi_blocks_u.push_back(std::make_pair(curr_img_pos,curr_block_len));
						clen -= curr_block_len;
						curr_img_pos += curr_block_len;
					}
				}
			}
		}
		/* phase 2) Break move-r blocks to avoid the images to overlap 
		 *          more than d consecutive blocks */
		break_phi_blocks(phi_blocks_u,M);
		/* phase 3) Fill in the move-r data structure blocks with block pointer, 
		            block offset, block starting point in the superblock, and flag
		            marking if the block is the last one in the superblock */
		r = phi_blocks_u.size() + phi_blocks_u.size()/M + M;
		block_w = 64 - __builtin_clzll(r);
		if(bit_width < (block_w + 2*uni_w + 1))
		{
			std::cerr << "	### Error! The current block vector of width " <<
			             bit_width << " cannot contain the move blocks" <<
			             " correctly. Please select a smaller superblock length or" <<
			             " let this function compute it automatically..." << std::endl;
			             exit(1);
		}
		else if(verbose)
		{		
			std::cout << "Memory displacement: |" << int(block_w) << " bits|" <<
			             int(uni_w) << " bits|" << int(uni_w) << " bits|1 bits|" <<
			             " total per block = " << int(block_w + 2*uni_w + 1) <<
			             " bits" << std::endl;
		}
		/* temp bitvector storing the block borders */
		// std::vector<usafe_t> superblock_c;
		{
			// bitvector borders;
			/* Compute the Elias-Fano predecessor data structure storing the block borders */
			{
				uint64_t val = 0, sum = 0;
				std::vector<uint64_t> b_borders{0}; b_borders.reserve(phi_blocks_u.size());
				for(const auto& p : phi_blocks_u)
				{
					sum += p.second;
					b_borders.push_back(sum);
				}
				borders.build(b_borders,sum);
			}
			uint64_t sum = 0, val = 0, sup = 0, b_count = 0;
			bool last;
			blocks.reserve(phi_blocks_u.size() + (phi_blocks_u.size()+(M-1))/M + M);
			// superblock_c.reserve(U/L + 1); superblock_c.push_back(sum);
			for(const auto& p : phi_blocks_u)
			{
				val = 0;
				b_count++;
				last = ((sum + p.second) % L) == 0;
				if(last) {
					sup++;
					//superblock_c.push_back(b_count);
				}

				auto b_pointer = borders.rank1(p.first+1)-1;
				auto next_block = borders.select1(b_pointer);
				val = ((((((
					  val | b_pointer) << uni_w) |
				      (p.first - next_block)) << uni_w) |
				      (sum % L)) << 1) |
				      static_cast<uint64_t>(last);

				blocks.push_back(val);
				sum += p.second;

				if(b_count % M == 0)
				{
					val = ((((((
						  static_cast<uint64_t>(0) | sup) << uni_w) |
					      static_cast<uint64_t>(0)) << uni_w) |
					      (sum % L)) << 1) |
					      static_cast<uint64_t>(0);
					blocks.push_back(val);
				}
			}
			// add block padding
			// superblock_c.push_back(b_count);
			int_t to_pad = M - (phi_blocks_u.size() % M);
			if(to_pad > 0)
			{
				for(size_t i = 0; i<to_pad; ++i){ blocks.push_back(0); }
				val = ((((((
					  static_cast<uint64_t>(0) | sup) << uni_w) |
				      static_cast<uint64_t>(0)) << uni_w) |
				      (sum % L)) << 1) |
				      static_cast<uint64_t>(0);
				blocks.push_back(val);
			}
		}
		/* print some stats */
		N = blocks.size();
		std::cout << "		- Phi-permutation move blocks = " << N << std::endl;
		std::cout << "		- Superblock length = " << L << std::endl;
		std::cout << "		- Universe size = " << U << std::endl;
	}

	inline uint64_t move(safe_t& p, safe_t& o, safe_t& c) const
	{
		/* compute correct position in the block vector
		   skipping the extra dummy blocks storing the 
		   superblock ids */
		safe_t b_p = p + p/M;
		/* compute the index of the dummy block containing
		   the next next superblock id */
		safe_t db_p = b_p + ((M+1) - ((b_p+1) % (M+1)));
		/* steps to walk in the permutation wrt the 
		   beginning of the pointed block */
		safe_t steps_to_walk = c + o;
		/* read one cache line in the buffer */
		__builtin_prefetch(&blocks[b_p], 0, 3);
		/* compute the vector containing the cumulative
		   counters of superblock shifts */
		safe_t shifts = blocks.at(b_p) & 1ULL;
		for(usafe_t i=1; i<(db_p-b_p)+1; ++i)
			{ shifts += (blocks.at(b_p+i) & 1ULL); }
		/* compute the resulting phi-query result */
		usafe_t res = L * ((blocks.at(db_p) >> (2*uni_w + 1)) - shifts) // superblock beginning
		              + ((blocks.at(b_p) >> 1) & ((1ULL << uni_w) - 1))  // offset wrt. superblock beginning
		              + o  // offset wrt. block beginning 
		              + c; // image offset
		/* walk forward in the permutation */
		safe_t cum_sum = 0, blen = 0, prev_blen = 0;
		safe_t prev_pos = steps_to_walk + 1;
		/* walk until the dummy node */
		for(usafe_t i=0; i<(db_p-b_p); ++i)
		{   /* compute cumulative block length */
			blen = (((blocks.at(b_p+i+1) >> 1) & ((1ULL << uni_w) - 1)) == 0) ?
			            L - ((blocks.at(b_p+i) >> 1) & ((1ULL << uni_w) - 1)) :
			            ((blocks.at(b_p+i+1) >> 1) & ((1ULL << uni_w) - 1)) -
			            ((blocks.at(b_p+i) >> 1) & ((1ULL << uni_w) - 1));
			cum_sum += blen;
			/* stop when we reach the correct blocks */ 
			if(((steps_to_walk+1) - cum_sum) < 0)
			{
				c = (prev_pos == 0  ? prev_blen : prev_pos) - 1;
				i = (prev_pos == 0) ? (i - 1) : i;
				/* update next pointer and offset */
				p = blocks.at(b_p+i) >> (uni_w * 2 + 1);
				o = (blocks.at(b_p+i) >> (uni_w + 1)) & ((1ULL << uni_w)-1);

				return res;
			}
			/* update previous positive residual positions count and block length */
			prev_pos = ((steps_to_walk+1) - cum_sum);
			prev_blen = blen;
		}
		/* walk until the last block */
		for(usafe_t i=(db_p-b_p)+1; i<M; ++i)
		{
			blen = (((blocks.at(b_p+i+1) >> 1) & ((1ULL << uni_w) - 1)) == 0) ?
			            L - ((blocks.at(b_p+i) >> 1) & ((1ULL << uni_w) - 1)) :
			            ((blocks.at(b_p+i+1) >> 1) & ((1ULL << uni_w) - 1)) -
			            ((blocks.at(b_p+i) >> 1) & ((1ULL << uni_w) - 1));
			cum_sum += blen;
			/* stop when we reach the correct blocks */ 
			if(((steps_to_walk+1) - cum_sum) < 0)
			{
				c = (prev_pos == 0 ? prev_blen : prev_pos) - 1;
				i = (prev_pos == 0)                   ? 
		            ((i == (db_p-b_p)+1) ? i-2 : i-1) :
							                       i;
				/* update next pointer and offset */
				p = blocks.at(b_p+i) >> (uni_w * 2 + 1);
				o = (blocks.at(b_p+i) >> (uni_w + 1)) & ((1ULL << uni_w)-1);

				return res;
			}
			/* update previous positive residual positions count and block length */
			prev_pos = ((steps_to_walk+1) - cum_sum);
			prev_blen = blen;
		}
		// process the last block
		c = (prev_pos == 0 ? prev_blen : prev_pos) - 1;
		safe_t i = (prev_pos == 0) ? ((M == (db_p-b_p)+1) ? M-2 : M-1) : M;
		/* update next pointer and offset */
		p = blocks.at(b_p+i) >> (uni_w * 2 + 1);
		o = (blocks.at(b_p+i) >> (uni_w + 1)) & ((1ULL << uni_w)-1);

		return res;
	}

	inline pair_t locate_block(safe_t i) const
	{
		safe_t block = borders.rank1(i+1)-1;
		safe_t offset = i - borders.select1(block);

		return std::make_pair(block+block/M,offset);
	}

	inline void get_next_block(safe_t& b, safe_t& o) const
	{
		uint64_t b_ = blocks.at(b);

		o = (b_ >> (1 + uni_w)) & ((1ULL << uni_w) - 1);
		b = (b_ >> (1 + 2*uni_w));
	}

	std::tuple<safe_t,safe_t,safe_t,bool> get_block(safe_t i)
	{
		usafe_t b = blocks.at(i);

		safe_t p  =  b >> (uni_w * 2 + 1);
		safe_t o  = (b >> (uni_w + 1)) & ((1ULL << uni_w)-1);
		safe_t c  = (b >> 1) & ((1ULL << uni_w)-1);
		bool   f  =  b & 1ULL;

		return make_tuple(p,o,c,f);
	}

	safe_t get_no_blocks() { return N; }

	safe_t get_d(){ return M; }

	safe_t get_universe_size(){ return U; }

	safe_t get_superblock_length(){ return L; }

	void load(std::istream& in)
	{
		in.read((char*)&uni_w, sizeof(uni_w));
		in.read((char*)&L, sizeof(L));
		in.read((char*)&N, sizeof(N));
		in.read((char*)&U, sizeof(U));

		blocks.resize(N);
		in.read(reinterpret_cast<char*>(blocks.data()), 
			    sizeof(typename bp_vector::value_type) * N);

		borders.load(in);
	}

	size_t serialize(std::ostream& out)
	{
		size_t w_bytes = 0;

		out.write((char*)&uni_w, sizeof(uni_w));
		out.write((char*)&L, sizeof(L));
		out.write((char*)&N, sizeof(N));
		out.write((char*)&U, sizeof(U));

		w_bytes += sizeof(uni_w) + sizeof(L) + sizeof(N) +
		           sizeof(U);

		w_bytes += sizeof(typename bp_vector::value_type) * N;
		out.write(reinterpret_cast<char*>(blocks.data()),
			      sizeof(typename bp_vector::value_type) * N);

		w_bytes += borders.serialize(out);

		return w_bytes;
	}
};
}

#endif 