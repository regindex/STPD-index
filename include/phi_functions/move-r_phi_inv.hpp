#ifndef MOVE_R_PHI_INV_HPP_
#define MOVE_R_PHI_INV_HPP_

//#include <elias_fano_sux.hpp>
//#include <sdsl/int_vector.hpp>
#include <move-r_large.hpp>
//#include <move-r.hpp>

namespace stpd{

class move_r_phi_inv
{
public:
	
	move_r_phi_inv(){} // empty constructor

	void build(const std::string bwt_filename, const std::string sa_filename, bool_t verbose = true)
	{
		std::ifstream bwt(bwt_filename,std::ifstream::binary);
		if(not bwt){ std::cerr << "Error opening the BWT file..." << std::endl; exit(1); }
		std::ifstream sa(sa_filename,std::ifstream::binary);
		if(not sa){ std::cerr << "Error opening the SA file..." << std::endl; exit(1); }

		bwt.seekg(0, bwt.end);
		usafe_t bwt_length = bwt.tellg();
		// skip first entry for $
		bwt.seekg(1, bwt.beg); sa.seekg(STORE_SIZE, sa.beg);

		std::vector<uint_t> phiPerm(bwt_length-1,0);

		usafe_t i=2,curr_sa=0,prev_sa=0,first_sa=0;
		//char curr=0; 
		//bwt.read(&curr,1);
		sa.read(reinterpret_cast<char*>(&first_sa), STORE_SIZE);
		prev_sa = first_sa;
		//std::cout << curr << " " << prev_sa-1 << std::endl;

		while(i++<bwt_length)
		{
			//bwt.read(&curr,1);
			sa.read(reinterpret_cast<char*>(&curr_sa), STORE_SIZE);

			//BWT[i] = curr;
			//SA[i++] = curr_sa-1;
			phiPerm[prev_sa-1] = curr_sa-1; 
			//std::cout << curr << " " << curr_sa-1 << std::endl;
			prev_sa = curr_sa;
		}

		//const uint_t MAX = (1ULL << sizeof(uint_t)*8) | (1ULL << sizeof(uint_t)*8) - 1;
		phiPerm[prev_sa-1] = first_sa-1;
		/*
		for(auto& e:phiPerm)
			std::cout << e << " ";
		std::cout << std::endl;
		*/
		using pair_t = std::pair<uint_t, uint_t>;
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
		/*
		for(auto& e:intervals)
			std::cout << "[" << e.first << "," << e.second << "] ";
		std::cout << std::endl;
		std::cout << "Number of runs = " << r << std::endl;
		std::cout << "Longest run = " << r_ << std::endl;
		*/
		bwt.close();
		sa.close();

		move.build(intervals,bwt_length-1,r,r_);
		/*
		init_phi2(13);
		usafe_t res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;
		res = phi_next2();
		std::cout << "######res = " << res << std::endl;

		exit(1);
		*/
	}

	void init_phi_(uint_t idx)
	{
		move.init_move(idx,qenv.b,qenv.o);
	}

	void init_phi(uint_t idx)
	{
		move.init_move2(idx,qenv2.b,qenv2.o,qenv2.o_);
	}

	int_t phi_next_()
	{ 
		move.move(qenv.b,qenv.o,qenv.c); // qenv.b: block id, qenv.o: offset, qenv.c: absolute pos

		return qenv.c;
	}

	int_t phi_next()
	{ 
		move.move2(qenv2.b,qenv2.o,qenv2.o_,qenv2.c); // qenv.b: block id, qenv.o: offset, qenv.c: absolute pos

		return qenv2.c;
	}

	int_t phi_unsafe(const uint_t idx) const
	{
		/*
		auto res = last.successor_rank(idx);

		return first[res.second-1] - (res.first - idx);
		*/
		return 0;
	}

	void load(std::istream& in)
	{
		//in.read((char*)&L, sizeof(L));
		move.load(in);
	}

	uint_t serialize(std::ostream& out)
	{
		uint_t w_bytes = 0;

		//out.write((char*)&L, sizeof(L));

		//w_bytes += sizeof(L);

		w_bytes += move.serialize(out);

		return w_bytes;
	}

private:
	class query_environment
	{
		public:
		usafe_t b; // block index
		usafe_t o; // block offset
		usafe_t c; // absolute position
	};

	class query_environment2
	{
		public:
		usafe_t b; // pointed block index
		usafe_t o; // pointed block offset
		usafe_t o_; // current block offset
		usafe_t c; // absolute position
	};
	
	move_r_large<> move;
	//move_r<> move;
	query_environment qenv;
	query_environment2 qenv2;
};
}

#endif // MOVE_PHI_INV_HPP_