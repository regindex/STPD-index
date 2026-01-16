add_library(text_oracle OBJECT
	bitpacked_text_oracle.hpp
	RLZ/RLZ_DNA.hpp
	RLZ/RLZ_DNA_sux.hpp
)
target_link_libraries(text_oracle PUBLIC common elias_fano sdsl divsufsort divsufsort64)