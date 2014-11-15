


#include "key.hpp"

void key_to_location(std::size_t key, integer* level, std::array<integer, NDIM>* loc) {
	for (integer d = 0; d != NDIM; ++d) {
		(*loc)[d] = 0;
	}
	*level = 0;
	while (key != 1) {
		for (integer d = NDIM-1; d>=0; --d) {
			(*loc)[d] <<= 1;
			(*loc)[d] |= (key & 1);
			key >>= 1;
		}
		(*level)++;
	}
}

std::size_t location_to_key(integer level, std::array<integer, NDIM> loc) {
	std::size_t key = 1;
	while (level > 0) {
		for (integer d = 0; d != NDIM; ++d) {
			key <<= 1;
			key |= (loc[d] & 1);
			loc[d] >>= 1;
		}
		level--;
	}
	return key;
}

