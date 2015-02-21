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

constexpr integer NX = 12;
constexpr integer FMM_BW = 2;
constexpr integer P = 0;
constexpr integer N3 = NX * NX * NX;
constexpr integer PP = P * P;
constexpr integer NDIM = 3;
constexpr integer NCHILD = 8;
constexpr integer NNEIGHBOR = 27;
constexpr integer MAXLEVEL = 5;





		;
constexpr integer NF = 2 * PP;

#endif /* DEF_HPP_ */
