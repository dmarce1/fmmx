/*
 * def.hpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#ifndef DEF_HPP_
#define DEF_HPP_

#include <cstdint>

using real = double;
using integer = std::int64_t;

constexpr integer NX = 8;
constexpr integer BW = 2;
constexpr integer P = 2;
constexpr integer N3 = NX * NX * NX;
constexpr integer PP = P * P;
constexpr integer NDIM = 3;
constexpr integer NCHILD = 8;
constexpr integer NNEIGHBOR = 27;
constexpr integer MAXLEVEL = 2;

#endif /* DEF_HPP_ */
