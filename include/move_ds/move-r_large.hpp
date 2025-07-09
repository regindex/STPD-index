#ifndef MOVE_R_LARGE_HPP_
#define MOVE_R_LARGE_HPP_

#include <elias_fano_sux.hpp>
#include <sdsl/int_vector.hpp>
#include "bitpacked_vector.hpp"

using pair_t = std::pair<uint_t, uint_t>;

/*
std::string uint128_to_string(__uint128_t value)
{
    if (value == 0) return "0";

    std::string result;
    while (value > 0) {
        result = "0123456789"[value % 10] + result;
        value /= 10;
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, __uint128_t value) {
    return os << uint128_to_string(value);
}
*/

namespace stpd{

template<class bitvector = sux::bits::EliasFano<>,class bp_vector = sdsl::int_vector<>>
class move_r_large
{
public:
	
	move_r_large(){} // empty constructor

	void build(const std::vector<pair_t>& intervals, usafe_t n, usafe_t r, usafe_t r_)
	{
		std::vector<uint64_t> borders(r,0);
		//block_begs.width(n_width);
		//block_begs.resize(r);

		usafe_t i, sum = 0;
		//block_begs[0] = sum;
		for(i=1;i<intervals.size();++i)
		{
			//std::cout << "i: " << i << std::endl;
			sum += intervals[i-1].second;
			//std::cout << "sum: " << sum << std::endl;
			borders[i] = sum;
			//block_begs[i] = sum;
		}
		//borders[r] = n;
		
		/* for(auto& b:borders)
			std::cout << b << " ";
		std::cout << std::endl; */

		block_borders.build(borders,n);

		n_width = 64 - __builtin_clzll(n);
		r_width = 64 - __builtin_clzll(r);
		r__width = 64 - __builtin_clzll(r_);
		n_width = 30;
		r_width = 25;
		r__width = 19;

		//std::cout << "widths = " << int(n_width) << "," << int(r_width) << "," << int(r__width) << std::endl;

		// 64 bits should be enough since phi runs are usually never larger than few milions positions
		if(n_width+r_width+r__width > 128)
		{
			std::cerr << "Blocks cannot be represented with 128 bits... exiting." << std::endl;
			exit(1);
		}
		blocks.width(n_width+r_width+r__width);
		blocks.resize_aligned(r+1);

		for(i=0;i<intervals.size();++i)
		{
			usafe_t rank   = block_borders.rank1(intervals[i].first+1)-1;
			usafe_t select = block_borders.select1(rank);
			__uint128_t bp_val = 0;
			bp_val = ((((bp_val | rank) << r__width) | (intervals[i].first - select)) << n_width) | borders[i];
			//std::cout << "interval (" << intervals[i].first << "," << intervals[i].first + intervals[i].second - 1 << ") -> ";
			//std::cout << rank << "," << intervals[i].first - select << "," << borders[i] << " - " << bp_val << std::endl;
			//blocks[i] = bp_val;
			blocks.set_aligned(i,bp_val);
			//break;
		}
		//blocks[r] = n;
		blocks.set_aligned(r,__uint128_t(n));

		/*
		//__uint128_t bp_val = blocks[0];
		__uint128_t bp_val = blocks.get_aligned(0);
		std::cout << "che esce= " << bp_val << std::endl;
		usafe_t c = bp_val & ((1ULL << n_width)-1);
		bp_val >>= n_width;
		usafe_t o = bp_val & ((1ULL << r__width)-1);
		usafe_t b = bp_val >> r__width;
		std::cout << b << " " << o << " " << c << std::endl;
		exit(1);
		*/

		return;
	}

	void init_move(usafe_t i,usafe_t& b,usafe_t& o)
	{
		std::cout << "widths = " <<  int(n_width) << " " << int(r_width) << " " << int(r__width) << std::endl;
		b = block_borders.rank1(i+1)-1;
		o = i - block_borders.select1(b);

		std::cout << "init-> " << b << "," << o << std::endl;
	}

	void init_move2(usafe_t i,usafe_t& b,usafe_t& o,usafe_t& o_)
	{
		b = block_borders.rank1(i+1)-1;
		o_ = i - block_borders.select1(b);

		//std::cout << "init-> " << b << "," << o << std::endl;

		__uint128_t bv_val = blocks.get_aligned(b); //__uint128_t bv_val = blocks[b];
		bv_val >>= n_width;
		o = bv_val & ((1ULL << r__width)-1);
		b = bv_val >> r__width;

		//std::cout << "init-> " << b << "," << o << "," << o_ << std::endl;
	}

	inline void move(usafe_t& b_,usafe_t& o_,usafe_t& c_)
	{
		
		std::cout << "		entra: " << b_ << " " << o_ << " " << c_ << std::endl;
		__uint128_t bv_val = blocks.get_aligned(b_); //__uint128_t bv_val = blocks[b_];
		std::cout << "bval= " << bv_val << std::endl;
		//usafe_t c = bv_val & ((1ULL << n_width)-1);
		bv_val >>= n_width;
		usafe_t o = bv_val & ((1ULL << r__width)-1);
		usafe_t b = bv_val >> r__width;
		
		std::cout << "[ " << b << " , " << o << std::endl;

		bv_val = blocks.get_aligned(b); //bv_val = blocks[b];
		usafe_t c = bv_val & ((1ULL << n_width)-1);
		c_ = c + o + o_;
		std::cout << "c_ = " << c_ << std::endl;

		__uint128_t next_bv_val = blocks.get_aligned(++b); // __uint128_t next_bv_val = blocks[++b];
		usafe_t next = next_bv_val & ((1ULL << n_width)-1);
		std::cout << "next = " << next << std::endl;
		while(next <= c_) 
		{
			c = next;
			next_bv_val = blocks.get_aligned(++b); //next_bv_val = blocks[++b];
			next = next_bv_val & ((1ULL << n_width)-1);
		}

		b_ = b - 1;
		o_ = c_ - c;
		std::cout << "-- " << b_ << " " << o_ << std::endl;
		std::cout << "risultato: " << c_ <<  std::endl;
	}

	inline void move2(usafe_t& b,usafe_t& o,usafe_t& o_,usafe_t& c)
	{
		__uint128_t bv_val = blocks.get_aligned(b); //__uint128_t bv_val = blocks[b];
		usafe_t c_ = bv_val & ((1ULL << n_width)-1);
		c = c_ + o + o_;
		////std::cout << "c = " << c << std::endl;

		__uint128_t next_bv_val = blocks.get_aligned(++b); // __uint128_t next_bv_val = blocks[++b];
		usafe_t next = next_bv_val & ((1ULL << n_width)-1);
		//std::cout << "next = " << next << std::endl;

		while(next <= c) 
		{
			c_ = next;
			bv_val = next_bv_val;
			next_bv_val = blocks.get_aligned(++b); //next_bv_val = blocks[++b];
			next = next_bv_val & ((1ULL << n_width)-1);
		}
		//b_ = b - 1;
		o_ = c - c_;
		bv_val >>= n_width;
		o = bv_val & ((1ULL << r__width)-1);
		b = bv_val >> r__width;

		//std::cout << b << "," << o << "," << o_ << std::endl;
	}

	/*
	inline usafe_t move_phi_next()
	{
		usafe_t bv_val = blocks[qenv.b];
		usafe_t b = bv_val >> r__width;
		usafe_t o = bv_val & ((1ULL << r__width)-1);

		//std::cout << "begin-> " << qenv.b << "," << qenv.o << std::endl;
		//std::cout << "block-> " << b << "," << o << std::endl;

		usafe_t current = block_borders.select1(b) + o;
		// Below it follows the branchless version of this cycle //
		while(qenv.o > 0)
		{
			usafe_t next = block_borders.select1(b+1);
			if((next - current) <= qenv.o)
			{
				qenv.o -= (next - current);
				o = 0;
				b++;
				current = next;
			}
			else
			{
				o += qenv.o;
				current += qenv.o;
				qenv.o = 0;
			}
		}
		*/
		/*
		while (qenv.o > 0)
		{
		    usafe_t next = block_borders.select1(b + 1);
		    usafe_t len = next - current;

		    // mask: 1 if len <= qenv.o, 0 otherwise
		    usafe_t cond = static_cast<usafe_t>(len <= qenv.o);

		    // conditional update
		    //usafe_t move = cond ? len : qenv.o;
		    usafe_t move = cond * len + (1 - cond) * qenv.o;

		    o       = (1 - cond) * (o + qenv.o);
		    current += move;
		    qenv.o  -= move;
		    b       += cond;
		}
		*/
		/*
		qenv.b = b;
		qenv.o = o;

		return current;
	}*/

	void load(std::istream& in)
	{
		in.read((char*)&n_width, sizeof(n_width));
		in.read((char*)&r_width, sizeof(r_width));
		in.read((char*)&r__width, sizeof(r__width));
		block_borders.load(in);
		blocks.load(in);
		//block_begs.load(in);
	}

	uint_t serialize(std::ostream& out)
	{
		uint_t w_bytes = 0;

		out.write((char*)&n_width, sizeof(n_width));
		out.write((char*)&r_width, sizeof(r_width));
		out.write((char*)&r__width, sizeof(r__width));

		w_bytes += sizeof(n_width) + sizeof(r_width) + sizeof(r__width);

		w_bytes += block_borders.serialize(out);
		w_bytes += blocks.serialize(out);
		//w_bytes += block_begs.serialize(out);

		return w_bytes;
	}

private:
	/*
	class query_environment
	{
		public:
		usafe_t b; // block index
		usafe_t o; // block offset
	};
	*/
	
	//move_r<> move;
	bitvector block_borders;
	bp_vector blocks_;
	bitpacked_vector blocks;
	//bp_vector block_begs;
	uint8_t n_width, r_width, r__width;
	//query_environment qenv;
};
}

#endif // MOVE_PHI_INV_HPP_