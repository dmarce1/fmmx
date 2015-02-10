/*
 * hydro.cpp
 *
 *  Created on: Feb 4, 2015
 *      Author: dmarce1
 */

#include "hydro.hpp"
#include "lane_emden.hpp"

constexpr real hydro_vars::ro_floor;

integer hydro_vars::ind3d(integer j, integer k, integer l, integer stride) {
	return l + stride * (k + stride * j);
}

void hydro_vars::store() {
	U0 = U;
}

real hydro_vars::cell_mass(integer i) const {
	return U[d0_i][i] * dx * dx * dx;
}

real hydro_vars::compute_du() {
	real amax = real(0);
	std::vector<std::vector<real>> F(nf_hydro, std::vector<real>(N3F, real(0)));
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer i = 0; i != N3; ++i) {
			dU[f][i] = real(0);
		}
	}
	std::vector<std::vector<real>> flx(nf_hydro, std::vector<real>(N3F, real(0)));
	for (integer d = 0; d != NDIM; ++d) {
		const integer dn[3] = { d == 0 ? 1 : 0, d == 1 ? 1 : 0, d == 2 ? 1 : 0 };
		for (integer j = 0; j != NX + dn[0]; ++j) {
			for (integer k = 0; k != NX + dn[1]; ++k) {
				for (integer l = 0; l != NX + dn[2]; ++l) {
					const auto vi = ind3d(j, k, l, NX + 1);
					const real rho_r = VR[d][ro_i][vi];
					const real rho_l = VL[d][ro_i][vi];
					const real v_r = VR[d][v0_i + d][vi];
					const real v_l = VL[d][v0_i + d][vi];
					const real vx_r = VR[d][vx_i][vi];
					const real vx_l = VL[d][vx_i][vi];
					const real vy_r = VR[d][vy_i][vi];
					const real vy_l = VL[d][vy_i][vi];
					const real vz_r = VR[d][vz_i][vi];
					const real vz_l = VL[d][vz_i][vi];
					const real pr_r = VR[d][pr_i][vi];
					const real pr_l = VL[d][pr_i][vi];
					const real et_r = pr_r / (fgamma - 1.0) + 0.5 * (vx_r * vx_r + vy_r * vy_r + vz_r * vz_r) * rho_r;
					const real et_l = pr_l / (fgamma - 1.0) + 0.5 * (vx_l * vx_l + vy_l * vy_l + vz_l * vz_l) * rho_l;
					const real a_r = std::sqrt(fgamma * std::max(pr_r, real(0)) / rho_r) + std::abs(v_r);
					const real a_l = std::sqrt(fgamma * std::max(pr_l, real(0)) / rho_l) + std::abs(v_l);
					const real a = std::max(a_r, a_l);
					amax = std::max(amax, a);
					flx[d0_i][vi] = real(0.5) * (rho_r * v_r + rho_l * v_l);
					flx[sx_i][vi] = real(0.5) * (rho_r * vx_r * v_r + rho_l * vx_l * v_l);
					flx[sy_i][vi] = real(0.5) * (rho_r * vy_r * v_r + rho_l * vy_l * v_l);
					flx[sz_i][vi] = real(0.5) * (rho_r * vz_r * v_r + rho_l * vz_l * v_l);
					flx[s0_i + d][vi] += real(0.5) * (pr_r + pr_l);
					flx[et_i][vi] = real(0.5) * ((et_r + pr_r) * v_r + (et_l + pr_l) * v_l);
					flx[d0_i][vi] -= real(0.5) * (rho_r - rho_l) * a;
					flx[sx_i][vi] -= real(0.5) * (vx_r * rho_r - vx_l * rho_l) * a;
					flx[sy_i][vi] -= real(0.5) * (vy_r * rho_r - vy_l * rho_l) * a;
					flx[sx_i][vi] -= real(0.5) * (vz_r * rho_r - vz_l * rho_l) * a;
					flx[et_i][vi] -= real(0.5) * (vz_r * et_r - vz_l * et_l) * a;
				}
			}
		}
		for (integer f = 0; f != nf_hydro; ++f) {
			for (integer j = 0; j != NX; ++j) {
				for (integer k = 0; k != NX; ++k) {
					for (integer l = 0; l != NX; ++l) {
						const auto iu0 = ind3d(j, k, l);
						const auto if0 = ind3d(j, k, l, NX + 1);
						const auto ifp = ind3d(j + dn[0], k + dn[1], l + dn[2], NX + 1);
						dU[f][iu0] -= (flx[f][ifp] - flx[f][if0]) / dx;
					}
				}
			}
		}
	}
	return cfl_factor * dx / amax;
}

hydro_vars::hydro_vars() {

	U = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
	dU = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
	U0 = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
	x = std::vector<real>(N3);
	y = std::vector<real>(N3);
	z = std::vector<real>(N3);
	r = std::vector<real>(N3);
	for (integer d = 0; d != NDIM; ++d) {
		VR[d] = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3F, real(0.0)));
		VL[d] = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3F, real(0.0)));
	}

}
void hydro_vars::dump_data(std::vector<double>::iterator& l, integer index) const {
	for (integer f = 0; f != nf_hydro; ++f) {
		*l++ = double(U[f][index]);
	}
}
void hydro_vars::initialize(real xcorner, real ycorner, real zcorner, real _dx) {
	dx = _dx;
	for (integer i = 0; i != N3; ++i) {
		auto const j = i / (NX * NX);
		auto const k = (i / NX) % NX;
		auto const l = i % NX;
		x[i] = xcorner + (real(j) + real(0.5)) * dx - real(0.5);
		y[i] = ycorner + (real(k) + real(0.5)) * dx - real(0.5);
		z[i] = zcorner + (real(l) + real(0.5)) * dx - real(0.5);
		r[i] = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);

	}
	/*
	 for (integer i = 0; i != N3; ++i) {
	 for (integer f = 0; f != nf_hydro; ++f) {
	 U[f][i] = 0.0;
	 }
	 real ro;
	 ro = lane_emden(8.0 * r[i]) / 8.0 / 8.0 / 8.0;
	 U[d0_i][i] = std::max(ro, ro_floor);

	 U[et_i][i] = real(4.0) * real(M_PI) / real(2.5) / (fgamma - real(1.0)) * std::pow(ro, real(5.0 / 3.0));

	 }
	 */
	for (integer i = 0; i != N3; ++i) {
		for (integer f = 0; f != nf_hydro; ++f) {
			U[f][i] = 0.0;
		}
		if (x[i] < 0.0) {
			U[d0_i][i] = 1.0;
			U[et_i][i] = 2.5;
		} else {
			U[d0_i][i] = 1.0e-1;
			U[et_i][i] = 0.125;
		}
	}
}

void hydro_vars::get_boundary(std::vector<real>::iterator i, integer d) const {
	constexpr integer xlb[2 * NDIM] = { 0, NX - bw, 0, 0, 0, 0 };
	constexpr integer ylb[2 * NDIM] = { 0, 0, 0, NX - bw, 0, 0 };
	constexpr integer zlb[2 * NDIM] = { 0, 0, 0, 0, 0, NX - bw };
	constexpr integer xub[2 * NDIM] = { bw, NX, NX, NX, NX, NX };
	constexpr integer yub[2 * NDIM] = { NX, NX, bw, NX, NX, NX };
	constexpr integer zub[2 * NDIM] = { NX, NX, NX, NX, bw, NX };
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer j = xlb[d]; j != xub[d]; ++j) {
			for (integer k = ylb[d]; k != yub[d]; ++k) {
				for (integer l = zlb[d]; l != zub[d]; ++l) {
					*i = U[f][ind3d(j, k, l)];
					++i;
				}
			}
		}
	}

}

bool hydro_vars::needs_refinement() const {
	return true;
}

void hydro_vars::update(real dt, integer rk, const std::vector<real>& phi, const std::vector<real>& gx,
		const std::vector<real>& gy, const std::vector<real>& gz) {
	real beta[2] = { 1.0, 0.5 };
	if (rk == 0) {
		U0 = U;
	}
	if (P != 0) {
		for (integer i = 0; i != N3; ++i) {
			dU[sx_i][i] += gx[i] * U[d0_i][i];
			dU[sy_i][i] += gy[i] * U[d0_i][i];
			dU[sz_i][i] += gz[i] * U[d0_i][i];
			dU[et_i][i] += gx[i] * U[sx_i][i];
			dU[et_i][i] += gy[i] * U[sy_i][i];
			dU[et_i][i] += gz[i] * U[sz_i][i];
		}
	}
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer i = 0; i != N3; ++i) {
			real u1 = U[f][i] + dU[f][i] * dt;
			real u0 = U0[f][i];
			U[f][i] = u1 * beta[rk] + u0 * (1.0 - beta[rk]);
		}
	}
	for (integer i = 0; i != N3; ++i) {
		U[d0_i][i] = std::max(U[d0_i][i], ro_floor);
	}
}

void hydro_vars::set_boundary(integer d, std::vector<real>::iterator* iter) {
	constexpr integer xoff[2 * NDIM] = { -bw, NX / 2 - bw, 0, 0, 0, 0 };
	constexpr integer yoff[2 * NDIM] = { 0, 0, -bw, NX / 2 - bw, 0, 0 };
	constexpr integer zoff[2 * NDIM] = { 0, 0, 0, 0, -bw, NX / 2 - bw };
	constexpr integer xbdim[2 * NDIM] = { NX / 2 + 2 * bw, NX / 2 + 2 * bw, NX, NX, NX, NX };
	constexpr integer ybdim[2 * NDIM] = { NX, NX, NX / 2 + 2 * bw, NX / 2 + 2 * bw, NX, NX };
	constexpr integer zbdim[2 * NDIM] = { NX, NX, NX, NX, NX / 2 + 2 * bw, NX / 2 + 2 * bw };
	constexpr integer ixlb[2 * NDIM] = { bw, 0, 0, 0, 0, 0 };
	constexpr integer ixub[2 * NDIM] = { NX / 2 + 2 * bw, NX / 2 + bw, NX, NX, NX, NX };
	constexpr integer iylb[2 * NDIM] = { 0, 0, bw, 0, 0, 0 };
	constexpr integer izlb[2 * NDIM] = { 0, 0, 0, 0, bw, 0 };
	constexpr integer iyub[2 * NDIM] = { NX, NX, NX / 2 + 2 * bw, NX / 2 + bw, NX, NX };
	constexpr integer izub[2 * NDIM] = { NX, NX, NX, NX, NX / 2 + 2 * bw, NX / 2 + bw };
	constexpr integer bxlb[2 * NDIM] = { 0, NX / 2 + bw, 0, 0, 0, 0 };
	constexpr integer bxub[2 * NDIM] = { bw, NX / 2 + 2 * bw, NX, NX, NX, NX };
	constexpr integer bylb[2 * NDIM] = { 0, 0, 0, NX / 2 + bw, 0, 0 };
	constexpr integer bzlb[2 * NDIM] = { 0, 0, 0, 0, 0, NX / 2 + bw };
	constexpr integer byub[2 * NDIM] = { NX, NX, bw, NX / 2 + 2 * bw, NX, NX };
	constexpr integer bzub[2 * NDIM] = { NX, NX, NX, NX, bw, NX / 2 + 2 * bw };
	constexpr auto this_n3 = (NX / 2 + 2 * bw) * NX * NX;
	constexpr auto slope_n3 = (NX + 2) * (NX + 2) * (NX + 2);
	std::vector<std::vector<real>> V(nf_hydro, std::vector<real>(this_n3, real(0)));
	std::vector<std::vector<real>> this_U(nf_hydro, std::vector<real>(this_n3, real(0)));
	const integer dn[3] = { d / 2 == 0 ? 1 : 0, d / 2 == 1 ? 1 : 0, d / 2 == 2 ? 1 : 0 };
	auto this_ind3d = [&](integer j, integer k, integer l) {
		return l + zbdim[d] * (k + ybdim[d] * j);
	};
	for (integer f = 0; f != nf_hydro; ++f) {
		if (iter) {
			for (integer j = bxlb[d]; j != bxub[d]; ++j) {
				for (integer k = bylb[d]; k != byub[d]; ++k) {
					for (integer l = bzlb[d]; l != bzub[d]; ++l) {
						const auto i = this_ind3d(j, k, l);
						this_U[f][i] = *(*iter)++;

					}
				}
			}
		} else {
			for (integer j = bxlb[d]; j != bxub[d]; ++j) {
				for (integer k = bylb[d]; k != byub[d]; ++k) {
					for (integer l = bzlb[d]; l != bzub[d]; ++l) {
						const auto ib = this_ind3d(j, k, l);
						integer iu;
						if (d == 0) {
							iu = ind3d(0, k, l);
						} else if (d == 1) {
							iu = ind3d(NX - 1, k, l);
						} else if (d == 2) {
							iu = ind3d(j, 0, l);
						} else if (d == 3) {
							iu = ind3d(j, NX - 1, l);
						} else if (d == 4) {
							iu = ind3d(j, k, 0);
						} else /*if( d == 5 ) */{
							iu = ind3d(j, k, NX - 1);
						}
						this_U[f][ib] = U[f][iu];
					}
				}
			}
		}
		for (integer j = ixlb[d]; j != ixub[d]; ++j) {
			for (integer k = iylb[d]; k != iyub[d]; ++k) {
				for (integer l = izlb[d]; l != izub[d]; ++l) {
					const auto ui = ind3d(j + xoff[d], k + yoff[d], l + zoff[d]);
					const auto vi = this_ind3d(j, k, l);
					this_U[f][vi] = U[f][ui];
				}
			}
		}
	}
	for (integer i = 0; i != this_n3; ++i) {
		const auto ro = std::max(this_U[d0_i][i], ro_floor);
		const auto et = this_U[et_i][i];
		const auto sx = this_U[sx_i][i];
		const auto sy = this_U[sy_i][i];
		const auto sz = this_U[sz_i][i];
		V[ro_i][i] = ro;
		V[vx_i][i] = sx / ro;
		V[vy_i][i] = sy / ro;
		V[vz_i][i] = sz / ro;
		V[pr_i][i] = (fgamma - 1.0) * (et - real(0.5) * (sx * sx + sy * sy + sz * sz) / ro);
	}
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer j = dn[0]; j != xbdim[d] - 2 * dn[0]; ++j) {
			for (integer k = dn[1]; k != ybdim[d] - 2 * dn[1]; ++k) {
				for (integer l = dn[2]; l != zbdim[d] - 2 * dn[2]; ++l) {
					const auto ri = ind3d(j + dn[0] + xoff[d], k + dn[1] + yoff[d], l + dn[2] + zoff[d], NX + 1);
					const auto vi = this_ind3d(j, k, l);
					VL[d / 2][f][ri] = V[f][vi];
				}
			}
		}
		for (integer j = 2 * dn[0]; j != xbdim[d] - dn[0]; ++j) {
			for (integer k = 2 * dn[1]; k != ybdim[d] - dn[1]; ++k) {
				for (integer l = 2 * dn[2]; l != zbdim[d] - dn[2]; ++l) {
					const auto ri = ind3d(j + xoff[d], k + yoff[d], l + zoff[d], NX + 1);
					const auto vi = this_ind3d(j, k, l);
					VR[d / 2][f][ri] = V[f][vi];
				}
			}
		}
	}
}
