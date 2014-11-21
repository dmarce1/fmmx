/*
 * key.hpp
 *
 *  Created on: Nov 14, 2014
 *      Author: dmarce1
 */

#ifndef KEY_HPP_
#define KEY_HPP_

#include "defs.hpp"
#include <array>
#include <hpx/hpx_fwd.hpp>



std::size_t location_to_key(integer level, std::array<integer, NDIM> loc);

void key_to_location(std::size_t, integer* level, std::array<integer, NDIM>*);

hpx::id_type key_to_locality(std::size_t key);


#endif /* KEY_HPP_ */
