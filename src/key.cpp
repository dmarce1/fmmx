#include "key.hpp"

void key_to_location(std::size_t key, integer* level, std::array<integer, NDIM>* loc) {
	for (integer d = 0; d != NDIM; ++d) {
		(*loc)[d] = 0;
	}
	*level = 0;
	while (key != 1) {
		for (integer d = 0; d < NDIM; ++d) {
			(*loc)[d] <<= 1;
			(*loc)[d] |= (key & 1);
			key >>= 1;
		}
		(*level)++;
	}
//	printf("%lx %li %li %li %li\n", key, *level, (*loc)[0], (*loc)[1], (*loc)[2]);
}


hpx::id_type key_to_locality(std::size_t key) {
	static auto locality_list = hpx::find_all_localities();
	static auto locality_cnt = locality_list.size();
	static std::size_t chunk_size = 64;
	return locality_list[(key / chunk_size) % locality_cnt];

}

std::size_t location_to_key(integer level, std::array<integer, NDIM> loc) {
	std::size_t key = 1;
	while (level > 0) {
		for (integer d = 0; d != NDIM; ++d) {
			key <<= 1;
			key |= (loc[NDIM - 1 - d] & 1);
			loc[NDIM - 1 - d] >>= 1;
		}
		level--;
	}
//	printf("%lx\n", key);
	return key;
}

