// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  stpd-index: implementation of the Suffix Tree Path Decomposition index
 */

#ifndef STPD_INDEX_HPP_
#define STPD_INDEX_HPP_

#include <chrono>
#include <malloc_count.h>
// inverse phi functions
#include <r-index_phi_inv.hpp>
#include <stpd-index_phi_inv.hpp>
// bitpacked text oracle
#include <bitpacked_text_oracle.hpp>
// stpd-array binary search
#include <stpd_array_binary_search.hpp>

namespace stpd{

template<class STPDArray, class textOracle, class phiFunction>
class stpd_index{

private:
	// phi function data structure
	phiFunction phi;
	// text oracle data structure
	textOracle O;
	// stpd array search data structure
	STPDArray S;

public:
	// empty constructor
	stpd_index(){}

	// build suffixient index by indexing the supermaximal extensions
	void build_colex_m(const std::string &text_filepath)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		O.build(text_filepath);
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(text_filepath+".colex_m",&O); 
		std::cout << "Step 3) Constructing the phi function..." << std::endl;
	  	phi.build(text_filepath+".rbwt",text_filepath+".pa");
	  	//phi.test_phi(3);
	  	
	  	std::cout << "Done!" << std::endl;
	}

	// build suffixient index by indexing the supermaximal extensions
	void build_colex_m_v2(const std::string &text_filepath)
	{
		std::cout << "Constructing the STPD-index for " << text_filepath << std::endl;
		std::cout << "Step 1) Constructing the random-access text oracle..." << std::endl;
		O.build(text_filepath);
		std::cout << "Step 2) Constructing the STPD-array binary search data structure..." << std::endl;
		S.build(text_filepath+".colex_m",&O); 
		std::cout << "Step 3) Constructing the phi function..." << std::endl;
	  	phi.build(text_filepath+".rbwt",text_filepath+".pa",S.get_array_point(),S.get_PA_size());
	  	//phi.test_phi(3);
	  	
	  	std::cout << "Done!" << std::endl;
	}
	
	usafe_t store(const std::string &index_filepath)
	{
		std::ofstream out(index_filepath);

		//usafe_t w_bytes = 0;
		//w_bytes += O.size();
		//w_bytes += S.store(out);
		usafe_t phi_size = phi.serialize(out);
		std::cout << "Storing phi data structure, size = " << phi_size << " bytes" << std::endl;
		
		out.close();

		return phi_size;
	}
	/*
	void load(const std::string &oracle_filepath, 
		      const std::string &index_filepath)
	{
		O.load(oracle_filepath);

		std::ifstream in(index_filepath);
		S.load(in,&O);
		in.close();
	}
	*/
};

}

#endif