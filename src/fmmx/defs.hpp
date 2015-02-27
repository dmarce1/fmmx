/*
 * def.hpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#ifndef DEF_HPP_
#define DEF_HPP_

#include <cstdint>

#ifdef MINI_HPX
#include "../hpx/hpx_fwd.hpp"
#include "../hpx/hpx.hpp"
#else
#include <hpx/hpx_init.hpp>
#include <hpx/include/components.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/lcos/when_any.hpp>
#include <hpx/hpx_fwd.hpp>
#include <hpx/lcos/local/mutex.hpp>
#include <hpx/runtime/actions/component_action.hpp>
#include <hpx/runtime/components/server/managed_component_base.hpp>

#endif
#include <boost/serialization/utility.hpp>

using real = double;
using integer = std::int64_t;

constexpr real FGAMMA = 7.0/4.0;
constexpr integer NMOM = 10;
constexpr integer NX = 8;
constexpr integer FMM_BW = 2;
constexpr integer HYDRO_NX = NX + 2;
constexpr integer P = 5;
constexpr integer N3 = NX * NX * NX;
constexpr integer PP = P * P;
constexpr integer NDIM = 3;
constexpr integer NCHILD = 8;
constexpr integer NNEIGHBOR = 27;
constexpr integer MAXLEVEL = 4;
constexpr integer NF = 2 * PP;

constexpr integer P000 = 0;
constexpr integer P001 = 1;
constexpr integer P010 = 2;
constexpr integer P100 = 3;
constexpr integer P110 = 4;
constexpr integer P101 = 5;
constexpr integer P011 = 6;
constexpr integer P200 = 7;
constexpr integer P020 = 8;
constexpr integer P002 = 9;


static constexpr std::array<integer, NNEIGHBOR> dir_z = { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
		1, 1, 1, 1, 1, 1, 1 };

static constexpr std::array<integer, NNEIGHBOR> dir_y = { -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1,
		-1, 0, 0, 0, 1, 1, 1 };

static constexpr std::array<integer, NNEIGHBOR> dir_x = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1,
		-1, 0, 1, -1, 0, 1 };


#endif /* DEF_HPP_ */
