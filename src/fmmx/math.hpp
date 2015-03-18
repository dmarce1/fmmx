/*
 * math.hpp
 *
 *  Created on: Mar 7, 2015
 *      Author: dmarce1
 */

#ifndef MATH_HPP_
#define MATH_HPP_

#include "defs.hpp"

constexpr integer MAX_GAUSS = 4;

const real quadpt41 = std::sqrt(real(3) / real(7) - real(2) / real(7) * std::sqrt(real(6) / real(5)));
const real quadpt42 = std::sqrt(real(3) / real(7) + real(2) / real(7) * std::sqrt(real(6) / real(5)));
const real quadwt41 = (real(18) + std::sqrt(30)) / real(36);
const real quadwt42 = (real(18) - std::sqrt(30)) / real(36);

static const real gauss_points[MAX_GAUSS][MAX_GAUSS] =
	{
		{ real(0) },
		{ -real(1) / std::sqrt(real(3)), +real(1) / std::sqrt(real(3)) },
		{ -std::sqrt(real(3) / real(5)), real(0), +std::sqrt(real(3) / real(5)) },
		{ -quadpt42, -quadpt41, +quadpt41, +quadpt42 } };

static const real gauss_vol_lb[MAX_GAUSS][MAX_GAUSS] =
	{
		{ -real(1) },
		{ -real(1), real(0) },
		{ -real(1), -real(4) / real(9), +real(4) / real(9) },
		{ -real(1), -real(1)+quadwt42, real(0), real(1) - quadwt42 } };

static const real gauss_weights[MAX_GAUSS][MAX_GAUSS] =
	{
		{ real(2) },
		{ real(1), real(1) },
		{ real(5) / real(9), real(8) / real(9), real(5) / real(9) },
		{ quadwt42, quadwt41, quadwt41, quadwt42 } };

real LegendreP(integer n, real x);
real dLegendreP_dx(integer n, real x);

real quadrature_1d(integer n, const std::function<real(real)>&);
real quadrature_2d(integer n, const std::function<real(real, real)>&);
real quadrature_3d(integer n, const std::function<real(real, real, real)>&);

#endif /* MATH_HPP_ */
