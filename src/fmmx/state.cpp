/*
 * state.cpp
 *
 *  Created on: Mar 7, 2015
 *      Author: dmarce1
 */

#include <cassert>

#include "defs.hpp"
#include "state.hpp"

constexpr real fgamma = 7.0 / 5.0;
constexpr real eps_floor = 1.0e-8;

real pressure_eos(real d, real eps) {
	return (fgamma - real(1)) * d * std::max(eps, eps_floor);
}

state con_to_characteristic(const state& U, const state& dU, integer dir) {
	const integer i1 = sxi + dir;
	const integer i2 = (sxi + (dir == 0 ? 1 : 0));
	const integer i3 = (sxi + (dir == 2 ? 1 : 2));
	const real drho = dU[d0i];
	const real rho = U[d0i];
	const real u = U[i1] / U[d0i];
	const real v = U[i2] / U[d0i];
	const real w = U[i3] / U[d0i];
	const real p = pressure_eos(rho,
			(U[egi] - rho * (u * u + v * v + w * w) / real(2)) / rho);
	const real c = std::sqrt(fgamma * p / rho);
	const real du = (dU[i1] - u * drho) / rho;
	const real dv = (dU[i2] - v * drho) / rho;
	const real dw = (dU[i3] - w * drho) / rho;
	const real dp = (fgamma - real(1))
			* (dU[egi] - rho * (u * du + v * dv + w * dw)
					- drho * (u * u + v * v + w * w) / real(2));
	state dC(HYDRO_NF);
	dC[0] = drho - dp / (c * c);
	dC[1] = dv;
	dC[2] = dw;
	dC[3] = du + dp / (rho * c);
	dC[4] = du - dp / (rho * c);
	return dC;
}

state characteristic_to_con(const state& U, const state& dC, integer dir) {
	const integer i1 = sxi + dir;
	const integer i2 = sxi + (dir == 0 ? 1 : 0);
	const integer i3 = sxi + (dir == 2 ? 1 : 2);
	const real rho = U[d0i];
	const real u = U[i1] / U[d0i];
	const real v = U[i2] / U[d0i];
	const real w = U[i3] / U[d0i];
	const real p = pressure_eos(rho,
			(U[egi] - rho * (u * u + v * v + w * w) / real(2)) / rho);
	const real c = std::sqrt(fgamma * p / rho);
	const real h = (p + U[egi]) / rho;

	state dU(HYDRO_NF);

	dU[0] = dC[0];
	dU[i1] = dC[0] * u;
	dU[i2] = dC[0] * v;
	dU[i3] = dC[0] * w;
	dU[4] = dC[0] * (u * u + v * v + w * w) / real(2);

	dU[0] += real(0);
	dU[i1] += real(0);
	dU[i2] += rho * dC[1];
	dU[i3] += real(0);
	dU[4] += rho * v * dC[1];

	dU[0] += real(0);
	dU[i1] += real(0);
	dU[i2] += real(0);
	dU[i3] += rho * dC[2];
	dU[4] += rho * w * dC[2];

	dU[0] += rho / (real(2) * c) * dC[3];
	dU[i1] += rho / (real(2) * c) * dC[3] * (u + c);
	dU[i2] += rho * v / (real(2) * c) * dC[3];
	dU[i3] += rho * w / (real(2) * c) * dC[3];
	dU[4] += rho / (real(2) * c) * dC[3] * (h + u * c);

	dU[0] += -rho / (real(2) * c) * dC[4];
	dU[i1] += -rho / (real(2) * c) * dC[4] * (u - c);
	dU[i2] += -rho * v / (real(2) * c) * dC[4];
	dU[i3] += -rho * w / (real(2) * c) * dC[4];
	dU[4] += -rho / (real(2) * c) * dC[4] * (h - u * c);

	return dU;
}

real v_(const state& u, integer dir) {
	const real s = u[sxi + dir];
	const real d = u[d0i];
	return s / d;
}

real ek_(const state& u) {
	real e = real(0);
	for (integer d = 0; d != NDIM; ++d) {
		e += std::pow(u[sxi + d], 2);
	}
	return e / (real(2) * u[d0i]);
}

real ei_(const state& u) {
	return u[egi] - ek_(u);
}

real p_(const state& u) {
	return pressure_eos(u[d0i], ei_(u) / u[d0i]);
}

real h_(const state& u) {
	return (p_(u) + u[egi]) / u[d0i];
}

real cs_(const state& u) {
	const real p = p_(u);
	const real d = u[d0i];
	return std::sqrt(fgamma * p / d);
}

state flux(const state& u, integer dir) {
	const real v = v_(u, dir);
	std::valarray<real> f = v * u;
	const real p = p_(u);
	f[sxi + dir] += p;
	f[egi] += v * p;
	return f;
}

real spectral_radius(const state& u, integer dir) {
	return std::abs(v_(u, dir)) + cs_(u);
}

state prim_to_con(const state& v) {
	const real D = v[roi];
	const real Sx = D * v[vxi];
	const real Sy = D * v[vyi];
	const real Sz = D * v[vzi];
	const real ek = (Sx * Sx + Sy * Sy + Sz * Sz) / (D * real(2));
	const real ei = D * (v[h0i] - ek) / fgamma;
	const real Eg = ei + ek;
	state u(HYDRO_NF);
	u[d0i] = D;
	u[sxi] = Sx;
	u[syi] = Sy;
	u[szi] = Sz;
	u[egi] = Eg;
	return u;
}

state con_to_prim(const state& u) {
	const real rho = u[d0i];
	const real vx = u[sxi] / rho;
	const real vy = u[syi] / rho;
	const real vz = u[szi] / rho;
	const real ek = rho * (vx * vx + vy * vy + vz * vz) / real(2);
	const real p = (fgamma - real(1)) * (u[egi] - ek);
	const real h = (p + u[egi]) / rho;
	state v(HYDRO_NF);
	v[roi] = rho;
	v[vxi] = vx;
	v[vyi] = vy;
	v[vzi] = vz;
	v[h0i] = h;
	return v;
}

state roe(const state& UR, const state& UL, integer dir) {
	const integer i1 = sxi + dir;
	const integer i2 = (sxi + (dir == 0 ? 1 : 0));
	const integer i3 = (sxi + (dir == 2 ? 1 : 2));

	state V0(HYDRO_NF), U0(HYDRO_NF);
	auto VR = con_to_prim(UR);
	auto VL = con_to_prim(UL);
	const auto wl = std::sqrt(VL[d0i]);
	const auto wr = std::sqrt(VR[d0i]);
	V0[d0i] = wl * wr;
	V0[vxi] = (wr * VR[vxi] + wl * VL[vxi]) / (wr + wl);
	V0[vyi] = (wr * VR[vyi] + wl * VL[vyi]) / (wr + wl);
	V0[vzi] = (wr * VR[vzi] + wl * VL[vzi]) / (wr + wl);
	V0[h0i] = (wr * VR[h0i] + wl * VL[h0i]) / (wr + wl);
	U0 = prim_to_con(V0);
	const auto a = cs_(U0);

	const auto u = V0[i1];
	const auto v = V0[i2];
	const auto w = V0[i3];
	const auto h = V0[h0i];
	const auto rho = V0[d0i];
	std::valarray<std::valarray<real>> E(std::valarray<real>(HYDRO_NF),
			HYDRO_NF);
	state lambda(HYDRO_NF);
	state delta_x(HYDRO_NF);

	const auto tmp = rho / (real(2) * a);

	const auto delta_p = p_(UR) - p_(UL);
	const auto delta_rho = VR[d0i] - VL[d0i];
	const auto delta_u = VR[i1] - VL[i1];
	const auto delta_v = VR[i2] - VL[i2];
	const auto delta_w = VR[i3] - VL[i3];

	lambda[0] = u;
	lambda[1] = u;
	lambda[2] = u;
	lambda[3] = u + a;
	lambda[4] = u - a;

	E[0][0] = real(1);
	E[1][0] = u;
	E[2][0] = v;
	E[3][0] = w;
	E[4][0] = (u * u + v * v + w * w) / real(2);

	E[0][1] = real(0);
	E[1][1] = real(0);
	E[2][1] = rho;
	E[3][1] = real(0);
	E[4][1] = rho * v;

	E[0][2] = real(0);
	E[1][2] = real(0);
	E[2][2] = real(0);
	E[3][2] = rho;
	E[4][2] = rho * w;

	E[0][3] = +tmp;
	E[1][3] = +tmp * (u + a);
	E[2][3] = +tmp * v;
	E[3][3] = +tmp * w;
	E[4][3] = +tmp * (h + u * a);

	E[0][4] = -tmp;
	E[1][4] = -tmp * (u - a);
	E[2][4] = -tmp * v;
	E[3][4] = -tmp * w;
	E[4][4] = -tmp * (h - u * a);

	delta_x[0] = delta_rho - delta_p / (a * a);
	delta_x[1] = delta_v;
	delta_x[2] = delta_w;
	delta_x[3] = delta_u + delta_p / (rho * a);
	delta_x[4] = delta_u - delta_p / (rho * a);

	auto phi0 = [](real lambda, real delta) {
		if( std::abs(lambda) < delta) {
			return (lambda*lambda + delta*delta)/(real(2.0)*delta);
		} else {
			return std::abs(lambda);
		}
	};

	state f(HYDRO_NF);
	f = (flux(UR, dir) + flux(UL, dir)) / real(2);
	for (int i = 0; i != HYDRO_NF; ++i) {
		for (int j = 0; j != HYDRO_NF; ++j) {
			f[i] -= E[i][j] * phi0(lambda[j], delta_u) * delta_x[j] / real(2);
		}
	}
	return f;
}

state Riemann_flux(const state& ur, const state& ul, integer dir) {
//	return roe(ur, ul, dir);
	const real a = std::max(spectral_radius(ur, dir), spectral_radius(ul, dir));
	return (flux(ur, dir) + flux(ul, dir) - a * (ur - ul)) / real(2);
}

