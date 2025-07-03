set(PHI_SOURCES r-index_phi_inv_sux.hpp r-index_phi_inv_intlv.hpp move-r_phi_inv.hpp)

add_library(phi_functions OBJECT ${PHI_SOURCES})
target_link_libraries(phi_functions PUBLIC bitvectors elias_fano move_ds)
