set(MOVE_SOURCES move-r.hpp move-r_large.hpp)

add_library(move_ds OBJECT ${MOVE_SOURCES})
target_link_libraries(move_ds PUBLIC bitvectors)
