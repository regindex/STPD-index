######### common #########

  add_library(common OBJECT
    common/common.hpp
  )
  target_include_directories(common PUBLIC common)

######### bitvector #########

  add_library(bitvector OBJECT
    bitvector/elias_fano_bitvector.hpp
    bitvector/succinct_bitvector.hpp
  )
  target_include_directories(bitvector PUBLIC bitvector)
  target_link_libraries(bitvector PUBLIC common sdsl divsufsort divsufsort64) 

######### elias_fano #########

  add_library(elias_fano OBJECT
    elias_fano/elias_fano_sux.hpp
    elias_fano/elias_fano_intlv.hpp
    elias_fano/elias_fano_static_sorted_map.hpp
  )
  target_include_directories(elias_fano PUBLIC elias_fano)
  target_link_libraries(elias_fano PUBLIC sux)

######### move #########

  add_library(move OBJECT
    move/uni_part_move-r.hpp
  )
  target_include_directories(move PUBLIC move)
  target_link_libraries(move PUBLIC elias_fano)

######### phi_function #########

  add_library(phi OBJECT
    phi_function/r-index_phi_inv_sux.hpp
    phi_function/r-index_phi_inv_intlv.hpp
  )
  target_include_directories(phi PUBLIC phi_function)
  target_link_libraries(phi PUBLIC bitvector elias_fano)

  add_library(move_phi OBJECT
    phi_function/move-r_phi_inv.hpp
  )
  target_include_directories(move_phi PUBLIC phi_function)
  target_link_libraries(move_phi PUBLIC elias_fano move)

######### text_oracle #########

add_subdirectory(text_oracle)
target_include_directories(text_oracle PUBLIC text_oracle)

######### binary_search #########

  add_library(binary_search OBJECT
    binary_search/stpd_array_binary_search.hpp
    binary_search/stpd_array_binary_search_opt.hpp
    binary_search/tabulated_binary_search_DNA.hpp
  )
  target_include_directories(binary_search PUBLIC binary_search)
  target_link_libraries(binary_search PUBLIC text_oracle bitvector)

######### text_index #########
  add_library(text_index OBJECT
    text_index/colex-index.hpp
  )
  target_include_directories(text_index PUBLIC text_index)
  target_link_libraries(text_index PUBLIC elias_fano move text_oracle binary_search)