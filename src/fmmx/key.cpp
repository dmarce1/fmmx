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
	static auto list = hpx::find_all_localities();
	integer level;
	std::array<integer, NDIM> loc;
	key_to_location(key, &level, &loc);
	std::size_t cnt = (std::size_t(1) << (3 * level));
	std::size_t index = key ^ cnt;
	double chunk_size = double(cnt) / double(list.size());
	if (chunk_size > 0.0) {
		index = std::min(std::size_t(double(index) / chunk_size), list.size() - 1);
	} else {
		index = 0;
	}
	return list[index];

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

