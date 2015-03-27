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

constexpr integer NX = 4;

//#define SOD_SHOCK
#define BLAST_WAVE

constexpr integer HYDRO_NX = NX + 2;
constexpr integer HYDRO_N3 = HYDRO_NX * HYDRO_NX * HYDRO_NX;
constexpr integer HYDRO_RK = 3;
constexpr integer HYDRO_P = 3;
constexpr integer HYDRO_PPP = (HYDRO_P + 2) * (HYDRO_P + 1) * (HYDRO_P) / 6;
constexpr integer HYDRO_NF = 5;

constexpr integer SILO_SUB_NX = HYDRO_P;
constexpr integer SILO_NX = NX * SILO_SUB_NX;
constexpr integer SILO_N3 = NX * NX * NX * SILO_SUB_NX * SILO_SUB_NX * SILO_SUB_NX;

constexpr integer FMM_NX = NX;
constexpr integer FMM_BW = 2;
constexpr integer FMM_P = 0;
constexpr integer FMM_N3 = FMM_NX * FMM_NX * FMM_NX;
constexpr integer FMM_PP = FMM_P * FMM_P;

constexpr integer NDIM = 3;
constexpr integer NCHILD = 8;
constexpr integer NNEIGHBOR = 27;
constexpr integer MAXLEVEL = 3;

static constexpr std::array<integer, NNEIGHBOR> dir_z =
	{ -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

static constexpr std::array<integer, NNEIGHBOR> dir_y =
	{ -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1 };

static constexpr std::array<integer, NNEIGHBOR> dir_x =
	{ -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1 };

#endif /* DEF_HPP_ */
