#ifndef MOVE_R_HPP_
#define MOVE_R_HPP_

#include <elias_fano_sux.hpp>
#include <sdsl/int_vector.hpp>

using pair_t = std::pair<uint_t, uint_t>;

namespace stpd{

template<class bitvector = sux::bits::EliasFano<>,class bp_vector = sdsl::int_vector<>>
class move_r
{
public:
	
	move_r(){} // empty constructor

	void build(const std::vector<pair_t>& intervals, usafe_t n, usafe_t r, usafe_t r_)
	{
		std::vector<uint64_t> borders(r,0);

		usafe_t i, sum = 0;
		for(i=1;i<intervals.size();++i)
		{
			//std::cout << "i: " << i << std::endl;
			sum += intervals[i-1].second;
			//std::cout << "sum: " << sum << std::endl;
			borders[i] = sum;
		}
		//borders[r] = n;
		
		/* for(auto& b:borders)
			std::cout << b << " ";
		std::cout << std::endl; */

		block_borders.build(borders,n);

		r_width = 64 - __builtin_clzll(r);
		r__width = 64 - __builtin_clzll(r_);

		//std::cout << "widths = " << int(r_width) << "," << int(r__width) << std::endl;

		blocks.width(r_width+r__width);
		blocks.resize(r);

		for(i=0;i<intervals.size();++i)
		{
			usafe_t rank   = block_borders.rank1(intervals[i].first+1)-1;
			usafe_t select = block_borders.select1(rank);
			usafe_t bp_val = ((0ULL | rank) << r__width) | (intervals[i].first - select);
			//std::cout << "interval (" << intervals[i].first << "," << intervals[i].first + intervals[i].second - 1 << ") -> ";
			//std::cout << rank << "," << intervals[i].first - select << " - " << bp_val << std::endl;
			blocks[i] = bp_val;
		}

		return;
	}

	void init_move(usafe_t i,usafe_t& b,usafe_t& o)
	{
		b = block_borders.rank1(i+1)-1;
		o = i - block_borders.select1(b);

		//std::cout << "init-> " << qenv.b << "," << qenv.o << std::endl;
	}

	inline void move(usafe_t& b_,usafe_t& o_,usafe_t& c_)
	{
		usafe_t bv_val = blocks[b_];
		usafe_t b = bv_val >> r__width;
		usafe_t o = bv_val & ((1ULL << r__width)-1);

		usafe_t current = block_borders.select1(b) + o;
		/* Below it follows the branchless version of this cycle */
		while(o_ > 0)
		{
			usafe_t next = block_borders.select1(b+1);
			if((next - current) <= o_)
			{
				o_ -= (next - current);
				o = 0;
				b++;
				current = next;
			}
			else
			{
				o += o_;
				current += o_;
				o_ = 0;
			}
		}

		b_ = b;
		o_ = o;
		c_ = current;
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
		in.read((char*)&r_width, sizeof(r_width));
		in.read((char*)&r__width, sizeof(r__width));
		block_borders.load(in);
		blocks.load(in);
	}

	uint_t serialize(std::ostream& out)
	{
		uint_t w_bytes = 0;

		out.write((char*)&r_width, sizeof(r_width));
		out.write((char*)&r__width, sizeof(r__width));

		w_bytes += sizeof(r_width) + sizeof(r__width);

		w_bytes += block_borders.serialize(out);
		w_bytes += blocks.serialize(out);

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
	bp_vector blocks;
	uint8_t r_width;
	uint8_t r__width;
	//query_environment qenv;
};
}

#endif // MOVE_PHI_INV_HPP_