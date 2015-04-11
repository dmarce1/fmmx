/*
 * state.hpp
 *
 *  Created on: Mar 7, 2015
 *      Author: dmarce1
 */

#ifndef STATE_HPP_
#define STATE_HPP_


#include "defs.hpp"

using state = std::valarray<real>;


constexpr real rho_floor = 1.0e-6;


constexpr integer d0i = 0;
constexpr integer sxi = 1;
constexpr integer syi = 2;
constexpr integer szi = 3;
constexpr integer egi = 4;

constexpr integer roi = 0;
constexpr integer vxi = 1;
constexpr integer vyi = 2;
constexpr integer vzi = 3;
constexpr integer h0i = 4;

state flux( const state&, integer, real );
state prim_to_con( const state& );
state con_to_prim( const state& );
state Riemann_flux(const state& ur, const state& ul, integer dir, real phi_r, real phi_l);
real spectral_radius(const state& u, integer dir);
state con_to_characteristic(const state& U, const state& dU, integer dir);
state characteristic_to_con(const state& U, const state& dC, integer dir);

real cs_(const state& u);


#endif /* STATE_HPP_ */
