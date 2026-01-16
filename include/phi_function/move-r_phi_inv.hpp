#ifndef MOVE_R_PHI_INV_HPP_
#define MOVE_R_PHI_INV_HPP_

#include <uni_part_move-r.hpp>
#include <fstream>

#define STORE_SIZE 5

namespace stpd{

class move_r_phi_inv
{
public:
	
	move_r_phi_inv(){} // empty constructor

	move_r_phi_inv(const std::string& bwt_filename, const std::string& sa_filename, bool verbose = true)
	{ build(bwt_filename,sa_filename,verbose); } 

	void build(const std::string& bwt_filename, const std::string& sa_filename, bool verbose = true)
	{
		std::ifstream bwt(bwt_filename,std::ifstream::binary);
		if(not bwt){ std::cerr << "Error opening the BWT file..." << std::endl; exit(1); }
		std::ifstream sa(sa_filename,std::ifstream::binary);
		if(not sa) { std::cerr << "Error opening the SA file..."  << std::endl; exit(1); }

		bwt.seekg(0, bwt.end);
		usafe_t bwt_length = bwt.tellg();
		// skip first entry for $
		bwt.seekg(1, bwt.beg); sa.seekg(STORE_SIZE, sa.beg);

		std::vector<uint_t> phiPerm(bwt_length-1,0);

		usafe_t i=2,curr_sa=0,prev_sa=0,first_sa=0;
		sa.read(reinterpret_cast<char*>(&first_sa), STORE_SIZE);

		prev_sa = first_sa;
		while(i++<bwt_length)
		{
			sa.read(reinterpret_cast<char*>(&curr_sa), STORE_SIZE);

			phiPerm[prev_sa-1] = curr_sa-1; 
			prev_sa = curr_sa;
		}

		phiPerm[prev_sa-1] = first_sa-1;

		using pair_t = std::pair<safe_t, safe_t>;
		std::vector<pair_t> intervals;

		uint_t start = 0, startPhi = phiPerm[0];
		uint_t r = 1, r_ = 0;
		for(i=1;i<phiPerm.size();++i)
		{
			if(phiPerm[i] != phiPerm[i-1]+1)
			{
				intervals.push_back(std::make_pair(startPhi,i-start));
				r_ = std::max(r_,static_cast<uint_t>(i-start));
				r++;
				start = i;
				startPhi = phiPerm[i];
			}
		}
		intervals.push_back(std::make_pair(startPhi,i-start));

		bwt.close();
		sa.close();

		/* input: move intervals and universe length */
		move.build(intervals,bwt_length-1);
	}

	safe_t phi(safe_t i) const
	{
		safe_t b, o, o_;
		auto sb = move.locate_block(i);
		b  = sb.first;
		o_ = sb.second;
		move.get_next_block(b,o);

		return move.move(b,o,o_);
	}

	inline void init_phi(safe_t i) const
	{
		qenv.sa = i;
		auto sb = move.locate_block(qenv.sa);
		qenv.b  = sb.first;
		qenv.o_ = sb.second;
		move.get_next_block(qenv.b,qenv.o);
	}

	inline safe_t phi_next() const
	{ 
		qenv.sa = move.move(qenv.b,qenv.o,qenv.o_);

		return qenv.sa;
	}

	void load(std::istream& in)
	{
		move.load(in);
	}

	uint_t serialize(std::ostream& out)
	{
		uint_t w_bytes = 0;

		w_bytes += move.serialize(out);

		return w_bytes;
	}

	uint_t get_bwt_size()
	{
		return move.get_universe_size();
	}

private:

	struct query_environment
	{
		safe_t b;  // pointed block index
		safe_t o;  // pointed block offset
		safe_t o_; // current block offset
		safe_t sa; // current sa position
	};
	
	uni_part_move_r<> move;
	mutable query_environment qenv;
};
}

#endif // MOVE_PHI_INV_HPP_