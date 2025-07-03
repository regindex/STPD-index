set(MOVE_SOURCES move-r.hpp)

add_library(move_ds OBJECT ${MOVE_SOURCES})
target_link_libraries(move_ds PUBLIC bitvectors)
