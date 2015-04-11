/*
 * hydro.cpp
 *
 *  Created on: Mar 6, 2015
 *      Author: dmarce1
 */

#include "hydro.hpp"
#include "math.hpp"
#include "exafmm.hpp"
#include "rk.hpp"
#include "lane_emden.hpp"

const integer di = HYDRO_NX * HYDRO_NX;
const integer dj = HYDRO_NX;
const integer dk = 1;

static exafmm_kernel exafmm;

const integer NGF = std::max(HYDRO_P, integer(1));
const integer NGC = std::max(HYDRO_P, integer(1));

state lane_emden_star(real x, real y, real z);
state blast_wave(real x, real y, real z);
state sod_shock(real x, real y, real z);

std::array<
		std::array<
				std::array<std::array<std::array<std::array<std::array<std::array<std::array<real, HYDRO_P>, HYDRO_P>, HYDRO_P>, HYDRO_P>, HYDRO_P>, HYDRO_P>, 3>,
				3>, 3> rad_inv;

std::function<state(real, real, real)> init_func = lane_emden_star;

state lane_emden_star(real x, real y, real z) {
	state u(HYDRO_NF);
	u = real(0);
	x -= real(1) / real(2);
	y -= real(1) / real(2);
	z -= real(1) / real(2);
	//x -= 0.1;
//	y -= 0.05;
	real alpha = 0.1;
	const real r = std::sqrt(x * x + y * y + z * z) / alpha;
	real theta = 0.0;
	real theta_floor = 1.0e-10;
	if (r < 4) {
		theta = lane_emden(r, .05);
	} else if (r < 1.0e-3) {
		theta = real(1);
	}
	theta = std::max(theta, theta_floor);
	u[d0i] = std::pow(theta, 1.5);
	const auto c0 = real(4) * real(M_PI) * alpha * alpha / (real(5) / real(3));
	u[egi] = std::pow(theta, 2.5) * c0;
	return u;
}

state blast_wave(real x, real y, real z) {
	state u(HYDRO_NF);
	u = real(0);
	x -= 0.5;
	y -= 0.5;
	z -= 0.5;
	real r2 = x * x + y * y + z * z;
	u[d0i] = 1.0;
	u[egi] = 1.0e+6 * exp(-r2 / 1.0e-3);
	return u;
}

state sod_shock(real x, real y, real z) {
	state u(HYDRO_NF);
	u = real(0);
	if (z > 0.5) {
		u[d0i] = 1.0;
		u[egi] = 2.5;
	} else {
		u[d0i] = 0.125;
		u[egi] = 0.25;
	}
	//printf("%e %e %e %e\n", x, y, z, u[0]);
	return u;
}

bool hydro::refinement_needed() const {
	for (integer i = 0; i != NX + 1; ++i) {
		for (integer j = 0; j != NX + 1; ++j) {
			for (integer k = 0; k != NX + 1; ++k) {
				const real x = x0 + dx * real(i) - real(1) / real(2);
				const real y = y0 + dx * real(j) - real(1) / real(2);
				const real z = z0 + dx * real(k) - real(1) / real(2);
				if (x * x + y * y + z * z < 0.25 * 0.25) {
					return true;
				}
			}
		}
	}
	return false;
}

void hydro::initialize() {

	const auto& gpt = gauss_points[HYDRO_P - 1];
	const auto& gwt = gauss_weights[HYDRO_P - 1];

	for (integer ii = 0; ii != HYDRO_N3; ++ii) {
		for (integer pp = 0; pp != HYDRO_PPP; ++pp) {
			U[0][ii][pp] = real(0);
		}
	}

	for (integer i = 0; i != HYDRO_NX; ++i) {
		for (integer j = 0; j != HYDRO_NX; ++j) {
			for (integer k = 0; k != HYDRO_NX; ++k) {
				integer ii = gindex(i, j, k);
				real x, y, z;
				x = x0 + real(i - .5) * dx;
				y = y0 + real(j - .5) * dx;
				z = z0 + real(k - .5) * dx;
				for (integer gx = 0; gx != HYDRO_P; ++gx) {
					for (integer gy = 0; gy != HYDRO_P; ++gy) {
						for (integer gz = 0; gz != HYDRO_P; ++gz) {
							const real x1 = gpt[gx] * dx / 2.0 + x;
							const real y1 = gpt[gy] * dx / 2.0 + y;
							const real z1 = gpt[gz] * dx / 2.0 + z;
							auto u = init_func(x1, y1, z1);
							for (integer l = 0; l != HYDRO_P; ++l) {
								real px = LegendreP(l, gpt[gx]);
								for (integer m = 0; m != HYDRO_P - l; ++m) {
									real py = LegendreP(m, gpt[gy]);
									for (integer n = 0; n != HYDRO_P - l - m; ++n) {
										real pz = LegendreP(n, gpt[gz]);
										const integer p = pindex(l, m, n);
										U[0][ii][p] += (px * py * pz / real(8)) * gwt[gx] * gwt[gy] * gwt[gz] * real(2 * l + 1) * real(2 * m + 1)
												* real(2 * n + 1) * u;
									}
								}
							}
						}
					}

				}
			}
		}
	}

}

std::vector<real> hydro::output_data() const {
	const auto& gpt = gauss_points[HYDRO_P - 1];
	std::vector<real> data(SILO_N3 * (HYDRO_NF + 4));
	auto iter = std::begin(data);
	for (integer i = 0; i != SILO_NX; ++i) {
		for (integer j = 0; j != SILO_NX; ++j) {
			for (integer k = 0; k != SILO_NX; ++k) {
				const integer ic = i % HYDRO_P;
				const integer ip = i / HYDRO_P;
				const integer jc = j % HYDRO_P;
				const integer jp = j / HYDRO_P;
				const integer kc = k % HYDRO_P;
				const integer kp = k / HYDRO_P;
				const real x = ((real(ip) + 0.5) + gpt[ic] * 0.5) / real(NX);
				const real y = ((real(jp) + 0.5) + gpt[jc] * 0.5) / real(NX);
				const real z = ((real(kp) + 0.5) + gpt[kc] * 0.5) / real(NX);
				const auto u = get_U_at(x, y, z);
				const auto g = get_gforce_at(x, y, z);
				for (integer f = 0; f != HYDRO_NF; ++f) {
					*iter++ = u[f];
				}
				*iter++ = g[3];
				*iter++ = g[0];
				*iter++ = g[1];
				*iter++ = g[2];
			}
		}
	}
	return data;
}

state hydro::value_at(const std::valarray<state>& u, real x, real y, real z) {
	real px, py, pz;
	state v(HYDRO_NF);
	v = real(0);
	for (integer l = 0; l != HYDRO_P; ++l) {
		px = LegendreP(l, x);
		for (integer m = 0; m != HYDRO_P - l; ++m) {
			py = LegendreP(m, y);
			for (integer n = 0; n != HYDRO_P - l - m; ++n) {
				pz = LegendreP(n, z);
				const integer p = pindex(l, m, n);
				v += u[p] * px * py * pz;
			}
		}
	}
	return v;
}

real hydro::derivative_at(const std::valarray<real>& u, real x, real y, real z, integer dir) {
	real px, py, pz;
	real v = real(0);
	for (integer l = 0; l != HYDRO_P; ++l) {
		px = dir != 0 ? LegendreP(l, x) : dLegendreP_dx(l, x);
		for (integer m = 0; m != HYDRO_P - l; ++m) {
			py = dir != 1 ? LegendreP(m, y) : dLegendreP_dx(m, y);
			for (integer n = 0; n != HYDRO_P - l - m; ++n) {
				pz = dir != 2 ? LegendreP(n, z) : dLegendreP_dx(n, z);
				const integer p = pindex(l, m, n);
				v += u[p] * px * py * pz;
			}
		}
	}
//	printf( "----%e\n", u[1]);
	return v;
}

real hydro::value_at(const std::valarray<real>& u, real x, real y, real z) {
	real px, py, pz;
	real v;
	v = real(0);
	for (integer l = 0; l != HYDRO_P; ++l) {
		px = LegendreP(l, x);
		for (integer m = 0; m != HYDRO_P - l; ++m) {
			py = LegendreP(m, y);
			for (integer n = 0; n != HYDRO_P - l - m; ++n) {
				pz = LegendreP(n, z);
				const integer p = pindex(l, m, n);
				v += u[p] * px * py * pz;
			}
		}
	}
	return v;
}

std::valarray<std::valarray<real> > valarray_2d_alloc(integer a, integer b) {
	std::valarray<std::valarray<real> > U;
	U.resize(a);
	for (integer i = 0; i != a; ++i) {
		U[i].resize(b);
	}
	return U;
}

std::valarray<std::valarray<std::valarray<real>>>valarray_3d_alloc
( integer a, integer b, integer c ) {
	std::valarray<std::valarray<std::valarray<real>>> U;
	U.resize(a);
	for( integer i = 0; i != a; ++i) {
		U[i] = valarray_2d_alloc(b,c);
	}
	return U;
}

std::valarray<std::valarray<std::valarray<std::valarray<real>>> >valarray_4d_alloc
( integer a, integer b, integer c, integer d ) {
	std::valarray<std::valarray<std::valarray<std::valarray<real>>>> U;
	U.resize(a);
	for( integer i = 0; i != a; ++i) {
		U[i] = valarray_3d_alloc(b,c,d);
	}
	return U;
}

integer hydro::pindex(integer n, integer m, integer l) {
	const integer lmn = l + m + n;
	const integer mn = n + m;
	return (lmn + 2) * (lmn + 1) * lmn / 6 + (mn + 1) * mn / 2 + n;
}

integer hydro::gindex(integer i, integer j, integer k, integer d) {
	return k + d * (j + d * i);

}

real hydro::next_du(integer rk, const std::vector<std::vector<real>>& fflux) {
	real amax = real(0);
	const real dxinv = real(2) / dx;
	const auto& gfpt = gauss_points[NGF - 1];
	const auto& gfwt = gauss_weights[NGF - 1];
	const auto& gcpt = gauss_points[NGC - 1];
	const auto& gcwt = gauss_weights[NGC - 1];
	for (integer i = 0; i != HYDRO_N3; ++i) {
		for (integer pp = 0; pp != HYDRO_PPP; ++pp) {
			dU[rk][i][pp] = real(0);
		}
	}
	for (integer ii = 1; ii != HYDRO_NX; ++ii) {
		for (integer jj = 1; jj != HYDRO_NX; ++jj) {
			for (integer kk = 1; kk != HYDRO_NX; ++kk) {
				const integer i = gindex(ii, jj, kk);
				for (integer g1 = 0; g1 != NGF; ++g1) {
					for (integer g2 = 0; g2 != NGF; ++g2) {
						const real x1 = gfpt[g1];
						const real x2 = gfpt[g2];
						const state urx = value_at(U[rk][i], -real(1), x1, x2);
						const state ury = value_at(U[rk][i], x1, -real(1), x2);
						const state urz = value_at(U[rk][i], x1, x2, -real(1));
						const state ulx = value_at(U[rk][i - di], +real(1), x1, x2);
						const state uly = value_at(U[rk][i - dj], x1, +real(1), x2);
						const state ulz = value_at(U[rk][i - dk], x1, x2, +real(1));
						const auto g_rx = gforce_at(psi[rk][i], U[rk][i], -real(1), x1, x2);
						const auto g_ry = gforce_at(psi[rk][i], U[rk][i], x1, -real(1), x2);
						const auto g_rz = gforce_at(psi[rk][i], U[rk][i], x1, x2, -real(1));
						const auto g_lx = gforce_at(psi[rk][i - di], U[rk][i - di], +real(1), x1, x2);
						const auto g_ly = gforce_at(psi[rk][i - dj], U[rk][i - dj], x1, +real(1), x2);
						const auto g_lz = gforce_at(psi[rk][i - dk], U[rk][i - dk], x1, x2, +real(1));
						state fx = Riemann_flux(urx, ulx, 0, g_rx[3], g_lx[3]);
						state fy = Riemann_flux(ury, uly, 1, g_ry[3], g_ly[3]);
						state fz = Riemann_flux(urz, ulz, 2, g_rz[3], g_lz[3]);
						if ((is_child_amr[i] == is_child_amr[i - di]) && (is_compute_cell[i] || is_compute_cell[i - di])) {
							amax = std::max(amax, spectral_radius(urx, 0));
							amax = std::max(amax, spectral_radius(ulx, 0));
							for (integer l = 0; l != HYDRO_P; ++l) {
								for (integer m = 0; m != HYDRO_P - l; ++m) {
									for (integer n = 0; n != HYDRO_P - l - m; ++n) {
										const integer p = pindex(l, m, n);
										const state ix = fx * (LegendreP(m, x1) * LegendreP(n, x2) * gfwt[g1] * gfwt[g2]);
										const real lsgn = real(l % 2 == 0 ? +1 : -1);
										const real factor = real((2 * l + 1) * (2 * m + 1) * (2 * n + 1)) / real(8);
										dU[rk][i][p] += ix * (lsgn * factor) * dxinv;
										dU[rk][i - di][p] -= ix * (factor * dxinv);
									}
								}
							}
						}
						if ((is_child_amr[i] == is_child_amr[i - dj]) && (is_compute_cell[i] || is_compute_cell[i - dj])) {
							amax = std::max(amax, spectral_radius(ury, 1));
							amax = std::max(amax, spectral_radius(uly, 1));
							for (integer l = 0; l != HYDRO_P; ++l) {
								for (integer m = 0; m != HYDRO_P - l; ++m) {
									for (integer n = 0; n != HYDRO_P - l - m; ++n) {
										const integer p = pindex(l, m, n);
										const state iy = fy * (LegendreP(l, x1) * LegendreP(n, x2) * gfwt[g1] * gfwt[g2]);
										const real msgn = real(m % 2 == 0 ? +1 : -1);
										const real factor = real((2 * l + 1) * (2 * m + 1) * (2 * n + 1)) / real(8);
										dU[rk][i][p] += iy * (msgn * factor) * dxinv;
										dU[rk][i - dj][p] -= iy * (factor * dxinv);
									}
								}
							}
						}
						if ((is_child_amr[i] == is_child_amr[i - dk]) && (is_compute_cell[i] || is_compute_cell[i - dk])) {
							amax = std::max(amax, spectral_radius(urz, 2));
							amax = std::max(amax, spectral_radius(ulz, 2));
							for (integer l = 0; l != HYDRO_P; ++l) {
								for (integer m = 0; m != HYDRO_P - l; ++m) {
									for (integer n = 0; n != HYDRO_P - l - m; ++n) {
										const integer p = pindex(l, m, n);
										const state iz = fz * (LegendreP(l, x1) * LegendreP(m, x2) * gfwt[g1] * gfwt[g2]);
										const real nsgn = real(n % 2 == 0 ? +1 : -1);
										const real factor = real((2 * l + 1) * (2 * m + 1) * (2 * n + 1)) / real(8);
										dU[rk][i][p] += iz * (nsgn * factor) * dxinv;
										dU[rk][i - dk][p] -= iz * (factor * dxinv);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	for (integer ci = 0; ci != NCHILD; ++ci) {
		auto iter = fflux[ci].begin();
		const integer cx = (ci >> 0) & 1;
		const integer cy = (ci >> 1) & 1;
		const integer cz = (ci >> 2) & 1;
		const integer xlb = 1 + cx * (NX / 2);
		const integer ylb = 1 + cy * (NX / 2);
		const integer zlb = 1 + cz * (NX / 2);
		const integer zub = zlb + NX / 2;
		const integer xub = xlb + NX / 2;
		const integer yub = ylb + NX / 2;
		state f(HYDRO_NF);
		integer fi;
		fi = 0 + cx;
		if (is_child_amr_face[fi]) {
			//printf("%li %li\n", ci, fflux[ci].size());
			const integer sgn = integer(2 * (fi % 2) - 1);
			const integer ii = 1 + (fi % 2) * (HYDRO_NX - 3);
			for (integer jj = ylb; jj != yub; ++jj) {
				for (integer kk = zlb; kk != zub; ++kk) {
					const integer i = gindex(ii, jj, kk);
					std::valarray<real> tmp(real(0), HYDRO_PPP * HYDRO_NF);
					std::copy(iter, iter + HYDRO_NF * HYDRO_PPP, std::begin(tmp));
					for (integer l = 0; l != HYDRO_P; ++l) {
						for (integer m = 0; m != HYDRO_P - l; ++m) {
							for (integer n = 0; n != HYDRO_P - l - m; ++n) {
								const integer p = pindex(l, m, n);
								std::copy(std::begin(tmp) + p * HYDRO_NF, std::begin(tmp) + (p + 1) * HYDRO_NF, std::begin(f));
								if ((is_child_amr[i] != is_child_amr[i - sgn * di])) {
									dU[rk][i][p] += real(sgn) * f * dxinv;
								}
							}
						}
					}
					iter += HYDRO_NF * HYDRO_PPP;
				}
			}
		}
		fi = 2 + cy;
		if (is_child_amr_face[fi]) {
			const integer sgn = integer(2 * (fi % 2) - 1);
			const integer jj = 1 + (fi % 2) * (HYDRO_NX - 3);
			for (integer ii = xlb; ii != xub; ++ii) {
				for (integer kk = zlb; kk != zub; ++kk) {
					const integer i = gindex(ii, jj, kk);
					std::valarray<real> tmp(real(0), HYDRO_PPP * HYDRO_NF);
					std::copy(iter, iter + HYDRO_NF * HYDRO_PPP, std::begin(tmp));
					for (integer l = 0; l != HYDRO_P; ++l) {
						for (integer m = 0; m != HYDRO_P - l; ++m) {
							for (integer n = 0; n != HYDRO_P - l - m; ++n) {
								const integer p = pindex(l, m, n);
								std::copy(std::begin(tmp) + p * HYDRO_NF, std::begin(tmp) + (p + 1) * HYDRO_NF, std::begin(f));
								if ((is_child_amr[i] != is_child_amr[i - sgn * dj])) {
									dU[rk][i][p] += real(sgn) * f * dxinv;
								}
							}
						}
					}
					iter += HYDRO_NF * HYDRO_PPP;
				}
			}
		}

		fi = 4 + cz;
		if (is_child_amr_face[fi]) {
			const integer sgn = integer(2 * (fi % 2) - 1);
			const integer kk = 1 + (fi % 2) * (HYDRO_NX - 3);
			for (integer ii = xlb; ii != xub; ++ii) {
				for (integer jj = ylb; jj != yub; ++jj) {
					const integer i = gindex(ii, jj, kk);
					std::valarray<real> tmp(real(0), HYDRO_PPP * HYDRO_NF);
					std::copy(iter, iter + HYDRO_NF * HYDRO_PPP, std::begin(tmp));
					for (integer l = 0; l != HYDRO_P; ++l) {
						for (integer m = 0; m != HYDRO_P - l; ++m) {
							for (integer n = 0; n != HYDRO_P - l - m; ++n) {
								const integer p = pindex(l, m, n);
								std::copy(std::begin(tmp) + p * HYDRO_NF, std::begin(tmp) + (p + 1) * HYDRO_NF, std::begin(f));
								if ((is_child_amr[i] != is_child_amr[i - sgn * dk])) {
									dU[rk][i][p] += real(sgn) * f * dxinv;
								}
							}
						}
					}
					iter += HYDRO_NF * HYDRO_PPP;
				}
			}
		}
	}
	for (integer ii = 1; ii != HYDRO_NX - 1; ++ii) {
		for (integer jj = 1; jj != HYDRO_NX - 1; ++jj) {
			for (integer kk = 1; kk != HYDRO_NX - 1; ++kk) {
				const integer i = gindex(ii, jj, kk);
				if (is_compute_cell[i]) {
					for (integer gx = 0; gx != NGC; ++gx) {
						for (integer gy = 0; gy != NGC; ++gy) {
							for (integer gz = 0; gz != NGC; ++gz) {
								const real x = gcpt[gx];
								const real y = gcpt[gy];
								const real z = gcpt[gz];
								const state u = value_at(U[rk][i], x, y, z);
								const auto g = gforce_at(psi[rk][i], U[rk][i], x, y, z);
								const real phi = g[3];
								const real fg_x = g[0];
								const real fg_y = g[1];
								const real fg_z = g[2];
								state fx = flux(u, 0, phi);
								state fy = flux(u, 1, phi);
								state fz = flux(u, 2, phi);
								amax = std::max(amax, spectral_radius(u, 0));
								amax = std::max(amax, spectral_radius(u, 1));
								amax = std::max(amax, spectral_radius(u, 2));
								for (integer l = 0; l != HYDRO_P; ++l) {
									for (integer m = 0; m != HYDRO_P - l; ++m) {
										for (integer n = 0; n != HYDRO_P - l - m; ++n) {
											const integer p = pindex(l, m, n);
											const real px = LegendreP(l, x);
											const real py = LegendreP(m, y);
											const real pz = LegendreP(n, z);
											const real dpx_dx = dLegendreP_dx(l, x);
											const real dpy_dy = dLegendreP_dx(m, y);
											const real dpz_dz = dLegendreP_dx(n, z);
											const real wt = gcwt[gx] * gcwt[gy] * gcwt[gz];
											const real factor = real((2 * l + 1) * (2 * m + 1) * (2 * n + 1)) / real(8);
											dU[rk][i][p] += fx * dpx_dx * py * pz * wt * factor * dxinv;
											dU[rk][i][p] += fy * px * dpy_dy * pz * wt * factor * dxinv;
											dU[rk][i][p] += fz * px * py * dpz_dz * wt * factor * dxinv;

											real rho0 = u[d0i];
											if (p == 0) {
												rho0 -= rho_floor;
											}
											dU[rk][i][p][sxi] += rho0 * fg_x * px * py * pz * wt * factor;
											dU[rk][i][p][syi] += rho0 * fg_y * px * py * pz * wt * factor;
											dU[rk][i][p][szi] += rho0 * fg_z * px * py * pz * wt * factor;
#ifndef TOTAL_ENERGY
											dU[rk][i][p][egi] += u[sxi] * fg_x * px * py * pz * wt * factor;
											dU[rk][i][p][egi] += u[syi] * fg_y * px * py * pz * wt * factor;
											dU[rk][i][p][egi] += u[szi] * fg_z * px * py * pz * wt * factor;
#endif
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
#ifdef TOTAL_ENERGY
	for (integer ii = 1; ii != HYDRO_NX - 1; ++ii) {
		for (integer jj = 1; jj != HYDRO_NX - 1; ++jj) {
			for (integer kk = 1; kk != HYDRO_NX - 1; ++kk) {
				const integer i = gindex(ii, jj, kk);
				if (is_compute_cell[i]) {
					for (integer gx = 0; gx != NGC; ++gx) {
						for (integer gy = 0; gy != NGC; ++gy) {
							for (integer gz = 0; gz != NGC; ++gz) {
								const real wt = gcwt[gx] * gcwt[gy] * gcwt[gz];
								const real x = gcpt[gx];
								const real y = gcpt[gy];
								const real z = gcpt[gz];
								const state du = value_at(dU[rk][i], x, y, z);
								const state u = value_at(U[rk][i], x, y, z);
								real rho0 = u[d0i];
								const real phi = gforce_at(psi[rk][i], U[rk][i], x, y, z)[3];
								for (integer l = 0; l != HYDRO_P; ++l) {
									for (integer m = 0; m != HYDRO_P - l; ++m) {
										for (integer n = 0; n != HYDRO_P - l - m; ++n) {
											const integer p = pindex(l, m, n);
											const real px = LegendreP(l, x) * real(2 * l + 1) / real(2);
											const real py = LegendreP(m, y) * real(2 * m + 1) / real(2);
											const real pz = LegendreP(n, z) * real(2 * n + 1) / real(2);
											dU[rk][i][p][egi] -= du[d0i] * phi * px * py * pz * wt * (real(1) - rho_floor / rho0);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
#endif
	if (amax == real(0)) {
		return std::numeric_limits<real>::max();
	} else {
		return dx / amax;
	}
}

gforce_t hydro::self_gforce(const std::valarray<std::valarray<real>>& U, real x, real y, real z) const {
	const auto& gpt = gauss_points[HYDRO_P - 1];
	const auto& gwt = gauss_weights[HYDRO_P - 1];

	gforce_t g;
	std::fill(std::begin(g), std::end(g), real(0));
	for (integer gx = 0; gx < HYDRO_P; ++gx) {
		for (integer gy = 0; gy < HYDRO_P; ++gy) {
			for (integer gz = 0; gz < HYDRO_P; ++gz) {
				const real rho0 = value_at(U, gpt[gx], gpt[gy], gpt[gz])[d0i];
				const real mass = rho0 * dx * dx * dx * gwt[gx] * gwt[gy] * gwt[gz] / real(8);
				const real x0 = (gpt[gx] - x)*dx;
				const real y0 = (gpt[gy] - y)*dx;
				const real z0 = (gpt[gz] - z)*dx;
				const real r2 = x0 * x0 + y0 * y0 + z0 * z0;
				const real r = std::sqrt(r2);
				if (r > 1.0e-10) {
					const real r3 = r * r2;
					g[3] -= mass / r;
					g[0] += x0 * mass / r3;
					g[1] += y0 * mass / r3;
					g[2] += z0 * mass / r3;
				}
			}
		}
	}
	return g;
}

std::vector<real> hydro::flux_correct_pack(integer rk) const {
	const auto& gfpt = gauss_points[NGF - 1];
	const auto& gfwt = gauss_weights[NGF - 1];
	std::vector<real> a;
	integer cnt = 0;
	for (integer i = 0; i != 2 * NDIM; ++i) {
		if (is_amr_face[i]) {
			++cnt;
		}
	}
	a.resize(cnt * HYDRO_NF * HYDRO_PPP * NX * NX / 4);
	auto iter = a.begin();
	for (integer face = 0; face != 2; ++face) {
		if (is_amr_face[face]) {
			integer ii = (face % 2 == 0) ? 3 : HYDRO_NX - 3;
			for (integer jp = 1; jp != HYDRO_NX - 1; jp += 2) {
				for (integer kp = 1; kp != HYDRO_NX - 1; kp += 2) {
					std::valarray<real> tmp(real(0), HYDRO_PPP * HYDRO_NF);
					for (integer jj = jp; jj < jp + 2; ++jj) {
						for (integer kk = kp; kk < kp + 2; ++kk) {
							const integer i = gindex(ii, jj, kk);
							const integer c1 = (jj - 1) % 2;
							const integer c2 = (kk - 1) % 2;
							for (integer g1 = 0; g1 != NGF; ++g1) {
								for (integer g2 = 0; g2 != NGF; ++g2) {
									const real x1 = (real(2 * c1) + gfpt[g1] - real(1)) / real(2);
									const real x2 = (real(2 * c2) + gfpt[g2] - real(1)) / real(2);
									const state urx = value_at(U[rk][i], -real(1), gfpt[g1], gfpt[g2]);
									const state ulx = value_at(U[rk][i - di], +real(1), gfpt[g1], gfpt[g2]);
									const real phi_rx = gforce_at(psi[rk][i], U[rk][i], -real(1), gfpt[g1], gfpt[g2])[3];
									const real phi_lx = gforce_at(psi[rk][i - di], U[rk][i - di], +real(1), gfpt[g1], gfpt[g2])[3];
									state fx = Riemann_flux(urx, ulx, 0, phi_rx, phi_lx);
									for (integer l = 0; l != HYDRO_P; ++l) {
										for (integer m = 0; m != HYDRO_P - l; ++m) {
											for (integer n = 0; n != HYDRO_P - l - m; ++n) {
												const state ix = fx * (LegendreP(m, x1) * LegendreP(n, x2) * gfwt[g1] * gfwt[g2]);
												const real lsgn = (face % 2 == 0) ? real(1) : real(l % 2 == 0 ? +1 : -1);
												const real factor = real((2 * l + 1) * (2 * m + 1) * (2 * n + 1)) / real(8);
												const state f = ix * lsgn * factor;
												for (integer ff = 0; ff != HYDRO_NF; ++ff) {
													tmp[HYDRO_NF * pindex(l, m, n) + ff] += f[ff] / real(4);
												}
											}
										}
									}
								}
							}
						}
					}
					std::copy(std::begin(tmp), std::end(tmp), iter);
					iter += HYDRO_PPP * HYDRO_NF;
				}
			}
		}
	}
	for (integer face = 2; face != 4; ++face) {
		if (is_amr_face[face]) {
			integer jj = (face % 2 == 0) ? 3 : HYDRO_NX - 3;
			for (integer ip = 1; ip != HYDRO_NX - 1; ip += 2) {
				for (integer kp = 1; kp != HYDRO_NX - 1; kp += 2) {
					std::valarray<real> tmp(real(0), HYDRO_PPP * HYDRO_NF);
					for (integer ii = ip; ii < ip + 2; ++ii) {
						for (integer kk = kp; kk < kp + 2; ++kk) {
							const integer i = gindex(ii, jj, kk);
							const integer c1 = (ii - 1) % 2;
							const integer c2 = (kk - 1) % 2;
							for (integer g1 = 0; g1 != NGF; ++g1) {
								for (integer g2 = 0; g2 != NGF; ++g2) {
									const real x1 = (real(2 * c1) + gfpt[g1] - real(1)) / real(2);
									const real x2 = (real(2 * c2) + gfpt[g2] - real(1)) / real(2);
									const state ury = value_at(U[rk][i], gfpt[g1], -real(1), gfpt[g2]);
									const state uly = value_at(U[rk][i - dj], gfpt[g1], +real(1), gfpt[g2]);
									const real phi_ry = gforce_at(psi[rk][i], U[rk][i], gfpt[g1], -real(1), gfpt[g2])[3];
									const real phi_ly = gforce_at(psi[rk][i - dj], U[rk][i - dj], gfpt[g1], +real(1), gfpt[g2])[3];
									state fy = Riemann_flux(ury, uly, 1, phi_ry, phi_ly);
									for (integer l = 0; l != HYDRO_P; ++l) {
										for (integer m = 0; m != HYDRO_P - l; ++m) {
											for (integer n = 0; n != HYDRO_P - l - m; ++n) {
												const state iy = fy * (LegendreP(l, x1) * LegendreP(n, x2) * gfwt[g1] * gfwt[g2]);
												const real msgn = (face % 2 == 0) ? real(1) : real(m % 2 == 0 ? +1 : -1);
												const real factor = real((2 * l + 1) * (2 * m + 1) * (2 * n + 1)) / real(8);
												const state f = iy * msgn * factor;
												for (integer ff = 0; ff != HYDRO_NF; ++ff) {
													tmp[HYDRO_NF * pindex(l, m, n) + ff] += f[ff] / real(4);
												}
											}
										}
									}
								}
							}
						}
					}
					std::copy(std::begin(tmp), std::end(tmp), iter);
					iter += HYDRO_PPP * HYDRO_NF;
				}
			}
		}
	}
	for (integer face = 4; face != 6; ++face) {
		if (is_amr_face[face]) {
			integer kk = (face % 2 == 0) ? 3 : HYDRO_NX - 3;
			for (integer ip = 1; ip != HYDRO_NX - 1; ip += 2) {
				for (integer jp = 1; jp != HYDRO_NX - 1; jp += 2) {
					std::valarray<real> tmp(real(0), HYDRO_PPP * HYDRO_NF);
					for (integer ii = ip; ii < ip + 2; ++ii) {
						for (integer jj = jp; jj < jp + 2; ++jj) {
							const integer i = gindex(ii, jj, kk);
							const integer c1 = (ii - 1) % 2;
							const integer c2 = (jj - 1) % 2;
							for (integer g1 = 0; g1 != NGF; ++g1) {
								for (integer g2 = 0; g2 != NGF; ++g2) {
									const real x1 = (real(2 * c1) + gfpt[g1] - real(1)) / real(2);
									const real x2 = (real(2 * c2) + gfpt[g2] - real(1)) / real(2);
									const state urz = value_at(U[rk][i], gfpt[g1], gfpt[g2], -real(1));
									const state ulz = value_at(U[rk][i - dk], gfpt[g1], gfpt[g2], +real(1));
									const real phi_rz = gforce_at(psi[rk][i], U[rk][i], gfpt[g1], gfpt[g2], -real(1))[3];
									const real phi_lz = gforce_at(psi[rk][i - dk], U[rk][i - dk], gfpt[g1], gfpt[g2], +real(1))[3];
									state fz = Riemann_flux(urz, ulz, 2, phi_rz, phi_lz);
									for (integer l = 0; l != HYDRO_P; ++l) {
										for (integer m = 0; m != HYDRO_P - l; ++m) {
											for (integer n = 0; n != HYDRO_P - l - m; ++n) {
												const state iz = fz * (LegendreP(l, x1) * LegendreP(m, x2) * gfwt[g1] * gfwt[g2]);
												const real nsgn = (face % 2 == 0) ? real(1) : real(n % 2 == 0 ? +1 : -1);
												const real factor = real((2 * l + 1) * (2 * m + 1) * (2 * n + 1)) / real(8);
												const state f = iz * nsgn * factor;
												for (integer ff = 0; ff != HYDRO_NF; ++ff) {
													tmp[HYDRO_NF * pindex(l, m, n) + ff] += f[ff] / real(4);
												}
											}
										}
									}

								}
							}
						}
					}
					std::copy(std::begin(tmp), std::end(tmp), iter);
					iter += HYDRO_PPP * HYDRO_NF;
				}
			}
		}
	}
	return a;
}

void hydro::enforce_physical_boundaries(integer rk, integer dir) {
	const integer dim = dir / 2;
	const integer si = dim + sxi;
	const real sgn = real(2 * (dir % 2) - 1);
	for (integer j = 0; j != HYDRO_NX; ++j) {
		for (integer k = 0; k != HYDRO_NX; ++k) {
			integer iii;
			if (dir == 0) {
				iii = gindex(0, j, k);
			} else if (dir == 1) {
				iii = gindex(HYDRO_NX - 1, j, k);
			} else if (dir == 2) {
				iii = gindex(j, 0, k);
			} else if (dir == 3) {
				iii = gindex(j, HYDRO_NX - 1, k);
			} else if (dir == 4) {
				iii = gindex(j, k, 0);
			} else /*if (dir == 5)*/{
				iii = gindex(j, k, HYDRO_NX - 1);
			}
			for (integer ppp = 1; ppp < HYDRO_PPP; ++ppp) {
				U[rk][iii][ppp] = real(0);
			}
			for (integer l = 0; l != HYDRO_P; ++l) {
				for (integer m = 0; m != HYDRO_P - l; ++m) {
					for (integer n = 0; n != HYDRO_P - l - m; ++n) {
						if (!((l != 0 && dim == 0) || (m != 0 && dim == 1) || (n != 0 && dim == 2))) {
							const integer p = pindex(l, m, n);
							if (dir == 0) {
								U[rk][iii][p] = U[rk][gindex(1, j, k)][p];
							} else if (dir == 1) {
								U[rk][iii][p] = U[rk][gindex(HYDRO_NX - 2, j, k)][p];
							} else if (dir == 2) {
								U[rk][iii][p] = U[rk][gindex(j, 1, k)][p];
							} else if (dir == 3) {
								U[rk][iii][p] = U[rk][gindex(j, HYDRO_NX - 2, k)][p];
							} else if (dir == 4) {
								U[rk][iii][p] = U[rk][gindex(j, k, 1)][p];
							} else if (dir == 5) {
								U[rk][iii][p] = U[rk][gindex(j, k, HYDRO_NX - 2)][p];
							}
						}
					}
				}
			}
			if (U[rk][iii][0][si] * sgn < real(0)) {
				U[rk][iii][0][egi] -= std::pow(U[rk][iii][0][si], 2) / U[rk][iii][0][d0i];
				U[rk][iii][0][si] = real(0);
			}
		}
	}
}

void hydro::enforce_physical_phi_boundaries(integer rk, integer dir) {
	for (integer l = 0; l != HYDRO_P; ++l) {
		for (integer m = 0; m != HYDRO_P - l; ++m) {
			for (integer n = 0; n != HYDRO_P - l - m; ++n) {
				for (integer j = 0; j != HYDRO_NX; ++j) {
					for (integer k = 0; k != HYDRO_NX; ++k) {
						const integer p = pindex(l, m, n);
						if (dir == 0) {
							psi[rk][gindex(0, j, k)][p] = psi[rk][gindex(1, j, k)][p];
						} else if (dir == 1) {
							psi[rk][gindex(HYDRO_NX - 1, j, k)][p] = psi[rk][gindex(HYDRO_NX - 2, j, k)][p];
						} else if (dir == 2) {
							psi[rk][gindex(j, 0, k)][p] = psi[rk][gindex(j, 1, k)][p];
						} else if (dir == 3) {
							psi[rk][gindex(j, HYDRO_NX - 1, k)][p] = psi[rk][gindex(j, HYDRO_NX - 2, k)][p];
						} else if (dir == 4) {
							psi[rk][gindex(j, k, 0)][p] = psi[rk][gindex(j, k, 1)][p];
						} else if (dir == 5) {
							psi[rk][gindex(j, k, HYDRO_NX - 1)][p] = psi[rk][gindex(j, k, HYDRO_NX - 2)][p];
						}
					}
				}
			}
		}
	}
}

void hydro::next_u(integer rk, real dt) {
	const auto& alpha = alpha_rk[HYDRO_RK - 1];
	const auto& beta = beta_rk[HYDRO_RK - 1];
	for (integer i = 0; i != HYDRO_N3; ++i) {
		for (integer p = 0; p != HYDRO_PPP; ++p) {
			U[rk + 1][i][p] = real(0);
			for (integer k = 0; k < rk + 1; ++k) {
				U[rk + 1][i][p] += alpha[rk][k] * U[k][i][p] + beta[rk][k] * dU[k][i][p] * dt;
			}
		}
	}
	if (rk == HYDRO_RK - 1) {
		U[0] = U[rk + 1];
	}
}

state hydro::get_U_at(real x, real y, real z) const {
	constexpr real dx = real(1) / real(HYDRO_NX - 2);
	const integer i = integer(x / dx) + 1;
	const integer j = integer(y / dx) + 1;
	const integer k = integer(z / dx) + 1;
	x -= (i - 1) * dx;
	y -= (j - 1) * dx;
	z -= (k - 1) * dx;
	x = real(2) * x / dx - real(1);
	y = real(2) * y / dx - real(1);
	z = real(2) * z / dx - real(1);
	return value_at(U[0][gindex(i, j, k)], x, y, z);

}

gforce_t hydro::get_gforce_at(real x, real y, real z) const {
	constexpr real dx = real(1) / real(HYDRO_NX - 2);
	const integer i = integer(x / dx) + 1;
	const integer j = integer(y / dx) + 1;
	const integer k = integer(z / dx) + 1;
	x -= (i - 1) * dx;
	y -= (j - 1) * dx;
	z -= (k - 1) * dx;
	x = real(2) * x / dx - real(1);
	y = real(2) * y / dx - real(1);
	z = real(2) * z / dx - real(1);
	return gforce_at(psi[0][gindex(i, j, k)], U[0][gindex(i, j, k)], x, y, z);

}

void hydro::restrict_unpack(integer rk, const std::vector<real>& data, integer ci) {
	const integer xlb = 1 + (NX / 2) * ((ci >> 0) & 1);
	const integer ylb = 1 + (NX / 2) * ((ci >> 1) & 1);
	const integer zlb = 1 + (NX / 2) * ((ci >> 2) & 1);
	auto iter = data.begin();
	for (integer i = xlb; i != xlb + NX / 2; ++i) {
		for (integer j = ylb; j != ylb + NX / 2; ++j) {
			for (integer k = zlb; k != zlb + NX / 2; ++k) {
				const integer ii = gindex(i, j, k);
				if (!is_child_amr[ii]) {
					for (integer pp = 0; pp != HYDRO_PPP; ++pp) {
						std::copy(iter, iter + HYDRO_NF, std::begin(U[rk][ii][pp]));
						iter += HYDRO_NF;
					}
				}
			}
		}
	}
}

std::vector<real> hydro::restrict_pack(integer rk) const {
	const auto& gpt = gauss_points[HYDRO_P - 1];
	const auto& gwt = gauss_weights[HYDRO_P - 1];
	std::vector<real> a(HYDRO_PPP * HYDRO_N3 / NCHILD * HYDRO_NF, real(0));

	auto iter = a.begin();

	for (integer ip = 1; ip != HYDRO_NX / 2; ++ip) {
		for (integer jp = 1; jp != HYDRO_NX / 2; ++jp) {
			for (integer kp = 1; kp != HYDRO_NX / 2; ++kp) {
				std::valarray<state> u(state(real(0), HYDRO_NF), HYDRO_PPP);
				if (!is_amr[gindex(2 * ip, 2 * jp, 2 * kp)]) {
					for (integer ic = 0; ic != 2; ++ic) {
						for (integer jc = 0; jc != 2; ++jc) {
							for (integer kc = 0; kc != 2; ++kc) {
								const integer ii = gindex(2 * ip + ic - 1, 2 * jp + jc - 1, 2 * kp + kc - 1);
								for (integer gx = 0; gx != HYDRO_P; ++gx) {
									const real x = (real(2 * ic) + gpt[gx] - real(1)) / real(2);
									for (integer gy = 0; gy != HYDRO_P; ++gy) {
										const real y = (real(2 * jc) + gpt[gy] - real(1)) / real(2);
										for (integer gz = 0; gz != HYDRO_P; ++gz) {
											const real z = (real(2 * kc) + gpt[gz] - real(1)) / real(2);
											for (integer l = 0; l < HYDRO_P; ++l) {
												const real px = LegendreP(l, x) * real(2 * l + 1) / real(4) * gwt[gx];
												for (integer m = 0; m < HYDRO_P - l; ++m) {
													const real py = LegendreP(m, y) * real(2 * m + 1) / real(4) * gwt[gy];
													for (integer n = 0; n < HYDRO_P - l - m; ++n) {
														const real pz = LegendreP(n, z) * real(2 * n + 1) / real(4) * gwt[gz];
														const integer p = pindex(l, m, n);
														u[p] += (px * py * pz) * value_at(U[rk][ii], gpt[gx], gpt[gy], gpt[gz]);
													}
												}
											}
										}
									}
								}
							}
						}
					}
					for (integer p = 0; p != HYDRO_PPP; ++p) {
						std::copy(std::begin(u[p]), std::end(u[p]), iter);
						iter += HYDRO_NF;
					}
				}
			}
		}
	}
	return a;
}

std::vector<real> hydro::pack_boundary(integer rk, integer d) const {
	const integer xlb = dir_x[d] == +1 ? HYDRO_NX - 2 : 1;
	const integer ylb = dir_y[d] == +1 ? HYDRO_NX - 2 : 1;
	const integer zlb = dir_z[d] == +1 ? HYDRO_NX - 2 : 1;
	const integer xub = dir_x[d] == -1 ? 2 : HYDRO_NX - 1;
	const integer yub = dir_y[d] == -1 ? 2 : HYDRO_NX - 1;
	const integer zub = dir_z[d] == -1 ? 2 : HYDRO_NX - 1;
	std::vector<real> a(HYDRO_NF * HYDRO_PPP * (xub - xlb) * (yub - ylb) * (zub - zlb));
	auto iter = a.begin();
	for (integer i = xlb; i < xub; ++i) {
		for (integer j = ylb; j < yub; ++j) {
			for (integer k = zlb; k < zub; ++k) {
				const integer ii = gindex(i, j, k);
				for (integer pp = 0; pp != HYDRO_PPP; ++pp) {
					for (integer f = 0; f != HYDRO_NF; ++f) {
						*iter++ = U[rk][ii][pp][f];
					}
				}
			}
		}
	}
	return a;
}

std::vector<real> hydro::pack_phi_boundary(integer rk, integer d) const {
	const integer xlb = dir_x[d] == +1 ? HYDRO_NX - 2 : 1;
	const integer ylb = dir_y[d] == +1 ? HYDRO_NX - 2 : 1;
	const integer zlb = dir_z[d] == +1 ? HYDRO_NX - 2 : 1;
	const integer xub = dir_x[d] == -1 ? 2 : HYDRO_NX - 1;
	const integer yub = dir_y[d] == -1 ? 2 : HYDRO_NX - 1;
	const integer zub = dir_z[d] == -1 ? 2 : HYDRO_NX - 1;
	std::vector<real> a(FMM_PP * (xub - xlb) * (yub - ylb) * (zub - zlb));
	auto iter = a.begin();
	for (integer i = xlb; i < xub; ++i) {
		for (integer j = ylb; j < yub; ++j) {
			for (integer k = zlb; k < zub; ++k) {
				const integer ii = gindex(i, j, k);
				for (integer pp = 0; pp != FMM_PP; ++pp) {
					*iter++ = psi[rk][ii][pp];
				}
			}
		}
	}
	return a;
}

void hydro::unpack_boundary(integer rk, const std::vector<real>& a, integer d) {
	integer xlb = 1;
	integer ylb = 1;
	integer zlb = 1;
	integer xub = HYDRO_NX - 1;
	integer yub = HYDRO_NX - 1;
	integer zub = HYDRO_NX - 1;
	if (dir_x[d] == -1) {
		xlb = 0;
		xub = 1;
	} else if (dir_x[d] == +1) {
		xlb = HYDRO_NX - 1;
		xub = HYDRO_NX;
	}
	if (dir_y[d] == -1) {
		ylb = 0;
		yub = 1;
	} else if (dir_y[d] == +1) {
		ylb = HYDRO_NX - 1;
		yub = HYDRO_NX;
	}
	if (dir_z[d] == -1) {
		zlb = 0;
		zub = 1;
	} else if (dir_z[d] == +1) {
		zlb = HYDRO_NX - 1;
		zub = HYDRO_NX;
	}
	auto iter = a.begin();
	for (integer i = xlb; i < xub; ++i) {
		for (integer j = ylb; j < yub; ++j) {
			for (integer k = zlb; k < zub; ++k) {
				const integer ii = gindex(i, j, k);
				for (integer pp = 0; pp != HYDRO_PPP; ++pp) {
					for (integer f = 0; f != HYDRO_NF; ++f) {
						U[rk][ii][pp][f] = *iter++;
					}
				}
			}
		}
	}
}

void hydro::unpack_phi_boundary(integer rk, const std::vector<real>& a, integer d) {
	integer xlb = 1;
	integer ylb = 1;
	integer zlb = 1;
	integer xub = HYDRO_NX - 1;
	integer yub = HYDRO_NX - 1;
	integer zub = HYDRO_NX - 1;
	if (dir_x[d] == -1) {
		xlb = 0;
		xub = 1;
	} else if (dir_x[d] == +1) {
		xlb = HYDRO_NX - 1;
		xub = HYDRO_NX;
	}
	if (dir_y[d] == -1) {
		ylb = 0;
		yub = 1;
	} else if (dir_y[d] == +1) {
		ylb = HYDRO_NX - 1;
		yub = HYDRO_NX;
	}
	if (dir_z[d] == -1) {
		zlb = 0;
		zub = 1;
	} else if (dir_z[d] == +1) {
		zlb = HYDRO_NX - 1;
		zub = HYDRO_NX;
	}
	auto iter = a.begin();
	for (integer i = xlb; i < xub; ++i) {
		for (integer j = ylb; j < yub; ++j) {
			for (integer k = zlb; k < zub; ++k) {
				const integer ii = gindex(i, j, k);
				for (integer pp = 0; pp != FMM_PP; ++pp) {
					psi[rk][ii][pp] = *iter++;
				}
			}
		}
	}
}

std::vector<real> hydro::pack_amr_boundary(integer rk, integer d, integer ci) const {
	const integer cx = (ci >> 0) & 1;
	const integer cy = (ci >> 1) & 1;
	const integer cz = (ci >> 2) & 1;
	const integer xlb = (dir_x[d] == +1 ? HYDRO_NX - 3 : 1);
	const integer ylb = (dir_y[d] == +1 ? HYDRO_NX - 3 : 1);
	const integer zlb = (dir_z[d] == +1 ? HYDRO_NX - 3 : 1);
	const integer xub = (dir_x[d] == -1 ? 3 : HYDRO_NX - 1);
	const integer yub = (dir_y[d] == -1 ? 3 : HYDRO_NX - 1);
	const integer zub = (dir_z[d] == -1 ? 3 : HYDRO_NX - 1);
	std::vector<real> data(HYDRO_PPP * NCHILD * (xub - xlb) * (yub - ylb) * (zub - zlb));
	auto iter = data.begin();
	const auto& gpt = gauss_points[HYDRO_P - 1];
	const auto& gwt = gauss_weights[HYDRO_P - 1];
	for (integer ic = xlb; ic < xub; ++ic) {
		for (integer jc = ylb; jc < yub; ++jc) {
			for (integer kc = zlb; kc < zub; ++kc) {
				real x, y, z;
				const integer ip = (ic - 1) / 2 + 1 + cx * NX / 2;
				const integer jp = (jc - 1) / 2 + 1 + cy * NX / 2;
				const integer kp = (kc - 1) / 2 + 1 + cz * NX / 2;
				const integer ii = gindex(ip, jp, kp);
				state uc(HYDRO_NF);
				for (integer l = 0; l < HYDRO_P; ++l) {
					for (integer m = 0; m < HYDRO_P - l; ++m) {
						for (integer n = 0; n < HYDRO_P - l - m; ++n) {
							uc = real(0);
							for (integer gx = 0; gx != HYDRO_P; ++gx) {
								real px = LegendreP(l, gpt[gx]) * real(2 * l + 1) / real(2);
								for (integer gy = 0; gy != HYDRO_P; ++gy) {
									real py = LegendreP(m, gpt[gy]) * real(2 * m + 1) / real(2);
									for (integer gz = 0; gz != HYDRO_P; ++gz) {
										real pz = LegendreP(n, gpt[gz]) * real(2 * n + 1) / real(2);
										x = gpt[gx] / real(2) + real((ic % 2) ^ 1) - real(1) / real(2);
										y = gpt[gy] / real(2) + real((jc % 2) ^ 1) - real(1) / real(2);
										z = gpt[gz] / real(2) + real((kc % 2) ^ 1) - real(1) / real(2);
										auto v = value_at(U[rk][ii], x, y, z);
										uc += (gwt[gx] * gwt[gy] * gwt[gz]) * px * py * pz * v;
									}
								}
							}
							std::copy(std::begin(uc), std::end(uc), iter);
							iter += HYDRO_NF;
						}
					}
				}
			}
		}
	}
	return data;
}

void hydro::unpack_amr_boundary(integer rk, const std::vector<real>& a, integer d) {
	const integer xlb = (dir_x[d] == +1 ? HYDRO_NX - 3 : 1);
	const integer ylb = (dir_y[d] == +1 ? HYDRO_NX - 3 : 1);
	const integer zlb = (dir_z[d] == +1 ? HYDRO_NX - 3 : 1);
	const integer xub = (dir_x[d] == -1 ? 3 : HYDRO_NX - 1);
	const integer yub = (dir_y[d] == -1 ? 3 : HYDRO_NX - 1);
	const integer zub = (dir_z[d] == -1 ? 3 : HYDRO_NX - 1);
	auto iter = a.begin();
	for (integer i = xlb; i < xub; ++i) {
		for (integer j = ylb; j < yub; ++j) {
			for (integer k = zlb; k < zub; ++k) {
				const integer ii = gindex(i, j, k);
				for (integer l = 0; l < HYDRO_P; ++l) {
					for (integer m = 0; m < HYDRO_P - l; ++m) {
						for (integer n = 0; n < HYDRO_P - l - m; ++n) {
							const integer p = pindex(l, m, n);
							for (integer f = 0; f != HYDRO_NF; ++f) {
								U[rk][ii][p][f] = *iter++;
							}
						}
					}
				}
			}
		}
	}
}

std::valarray<real> minmod(real a, real b) {
	state c(HYDRO_NF);
	constexpr auto half = real(1) / real(2);
	c = (std::copysign(half, a) + std::copysign(half, b)) * std::min(std::abs(a), std::abs(b));
	return c;
}
;

std::valarray<real> minmaxmod(real a, real b1, real b2) {
	real c;
	for (integer i = 0; i != HYDRO_NF; ++i) {
		if (a > real(0)) {
			if (b1 > real(0) && b2 > real(0)) {
				c = std::max(b1, b2);
			} else if (b1 > real(0)) {
				c = b1;
			} else if (b2 > real(0)) {
				c = b2;
			} else {
				c = real(0);
			}
		} else {
			if (b1 < real(0) && b2 < real(0)) {
				c = std::min(b1, b2);
			} else if (b1 < real(0)) {
				c = b1;
			} else if (b2 < real(0)) {
				c = b2;
			} else {
				c = real(0);
			}
		}

	}
	return minmod(a, c);
}
;

std::valarray<real> minmod(const std::valarray<real>& a, const std::valarray<real>& b) {
	state c(HYDRO_NF);
	constexpr auto half = real(1) / real(2);
	for (integer i = 0; i != HYDRO_NF; ++i) {
		c[i] = (std::copysign(half, a[i]) + std::copysign(half, b[i])) * std::min(std::abs(a[i]), std::abs(b[i]));
	}
	return c;
}
;

std::valarray<real> minmaxmod(const std::valarray<real>& a, const std::valarray<real>& b1, const std::valarray<real> b2) {
	state c(HYDRO_NF);
	for (integer i = 0; i != HYDRO_NF; ++i) {
		if (a[i] > real(0)) {
			if (b1[i] > real(0) && b2[i] > real(0)) {
				c[i] = std::max(b1[i], b2[i]);
			} else if (b1[i] > real(0)) {
				c[i] = b1[i];
			} else if (b2[i] > real(0)) {
				c[i] = b2[i];
			} else {
				c[i] = real(0);
			}
		} else {
			if (b1[i] < real(0) && b2[i] < real(0)) {
				c[i] = std::min(b1[i], b2[i]);
			} else if (b1[i] < real(0)) {
				c[i] = b1[i];
			} else if (b2[i] < real(0)) {
				c[i] = b2[i];
			} else {
				c[i] = real(0);
			}
		}

	}
	return minmod(a, c);
}
;

bool real_eq(real a, real b) {
	constexpr real delta = 1.0e-11;
	if (a == real(0) && b == real(0)) {
		return true;
	}
	return std::abs(a - b) / (std::abs(a) + std::abs(b)) < delta;
}

void hydro::apply_limiter(integer rk) {

	auto& u = U[rk];
	auto ux = u;
	auto uy = u;
	auto uz = u;
	auto u_tilde = u;

	for (integer i = di; i < HYDRO_N3 - di; ++i) {
		for (integer m = 0; m < HYDRO_P - 1; ++m) {
			for (integer n = 0; n < HYDRO_P - m - 1; ++n) {
				for (integer l = HYDRO_P - 2 - n - m; l >= 0; --l) {
					const integer p = pindex(l, m, n);
					const integer pp1 = pindex(l + 1, m, n);
					state dup = (u[i + di][p] - u[i][p]) / real(2 * l + 1);
					state dum = (u[i][p] - u[i - di][p]) / real(2 * l + 1);
					state dur = (u[i + di][p] - u[i][p]) / real(2 * l + 1) - u[i + di][pp1];
					state dul = (u[i][p] - u[i - di][p]) / real(2 * l + 1) - u[i - di][pp1];
					state du0 = u[i][pp1];
					dup = con_to_characteristic(u[i][0], dup, 0);
					dum = con_to_characteristic(u[i][0], dum, 0);
					dur = con_to_characteristic(u[i][0], dur, 0);
					dul = con_to_characteristic(u[i][0], dul, 0);
					du0 = con_to_characteristic(u[i][0], du0, 0);
					state u1 = minmod(dup, dum);
					state u2 = minmod(dur, dul);
					u_tilde[i][pp1] = characteristic_to_con(u[i][0], minmaxmod(du0, u1, u2), 0);
				}
				for (integer f = 0; f != HYDRO_NF; ++f) {
					for (integer l = HYDRO_P - 2 - n - m; l >= 0; --l) {
						const integer pp1 = pindex(l + 1, m, n);
						if (!real_eq(u_tilde[i][pp1][f], u[i][pp1][f])) {
							ux[i][pp1][f] = u_tilde[i][pp1][f];
						} else {
							break;
						}
					}
				}
			}
		}
	}

	for (integer i = dj; i < HYDRO_N3 - dj; ++i) {
		for (integer l = 0; l < HYDRO_P - 1; ++l) {
			for (integer n = 0; n < HYDRO_P - l - 1; ++n) {
				for (integer m = HYDRO_P - 2 - l - n; m >= 0; --m) {
					const integer p = pindex(l, m, n);
					const integer pp1 = pindex(l, m + 1, n);
					state dup = (u[i + dj][p] - u[i][p]) / real(2 * m + 1);
					state dum = (u[i][p] - u[i - dj][p]) / real(2 * m + 1);
					state du0 = u[i][pp1];
					state dur = (u[i + dj][p] - u[i][p]) / real(2 * m + 1) - u[i + dj][pp1];
					state dul = (u[i][p] - u[i - dj][p]) / real(2 * m + 1) - u[i - dj][pp1];
					dup = con_to_characteristic(u[i][0], dup, 1);
					dum = con_to_characteristic(u[i][0], dum, 1);
					dur = con_to_characteristic(u[i][0], dur, 1);
					dul = con_to_characteristic(u[i][0], dul, 1);
					du0 = con_to_characteristic(u[i][0], du0, 1);
					state u1 = minmod(dup, dum);
					state u2 = minmod(dur, dul);
					u_tilde[i][pp1] = characteristic_to_con(u[i][0], minmaxmod(du0, u1, u2), 1);
				}
				for (integer f = 0; f != HYDRO_NF; ++f) {
					for (integer m = HYDRO_P - 2 - l - n; m >= 0; --m) {
						const integer pp1 = pindex(l, m + 1, n);
						if (!real_eq(u_tilde[i][pp1][f], u[i][pp1][f])) {
							uy[i][pp1][f] = u_tilde[i][pp1][f];
						} else {
							break;
						}
					}
				}
			}
		}
	}

	for (integer i = dk; i < HYDRO_N3 - dk; ++i) {
		for (integer m = 0; m < HYDRO_P - 1; ++m) {
			for (integer l = 0; l < HYDRO_P - m - 1; ++l) {
				for (integer n = HYDRO_P - 2 - m - l; n >= 0; --n) {
					const integer p = pindex(l, m, n);
					const integer pp1 = pindex(l, m, n + 1);
					state dup = (u[i + dk][p] - u[i][p]) / real(2 * n + 1);
					state dum = (u[i][p] - u[i - dk][p]) / real(2 * n + 1);
					state du0 = u[i][pp1];
					state dur = (u[i + dk][p] - u[i][p]) / real(2 * n + 1) - u[i + dk][pp1];
					state dul = (u[i][p] - u[i - dk][p]) / real(2 * n + 1) - u[i - dk][pp1];
					dup = con_to_characteristic(u[i][0], dup, 2);
					dum = con_to_characteristic(u[i][0], dum, 2);
					dur = con_to_characteristic(u[i][0], dur, 2);
					dul = con_to_characteristic(u[i][0], dul, 2);
					du0 = con_to_characteristic(u[i][0], du0, 2);
					state u1 = minmod(dup, dum);
					state u2 = minmod(dur, dul);
					u_tilde[i][pp1] = characteristic_to_con(u[i][0], minmaxmod(du0, u1, u2), 2);
				}
				for (integer f = 0; f != HYDRO_NF; ++f) {
					for (integer n = HYDRO_P - 2 - m - l; n >= 0; --n) {
						const integer pp1 = pindex(l, m, n + 1);
						if (!real_eq(u_tilde[i][pp1][f], u[i][pp1][f])) {
							uz[i][pp1][f] = u_tilde[i][pp1][f];
						} else {
							break;
						}
					}
				}
			}
		}
	}

	for (integer i = 0; i < HYDRO_N3; ++i) {
		for (integer pp = 1; pp != HYDRO_PPP; ++pp) {
			u[i][pp] = minmod(ux[i][pp], minmod(uz[i][pp], uy[i][pp]));
		}
	}

	for (integer ii = 0; ii != HYDRO_N3; ++ii) {
		const integer GP = NGC;
		const auto& gpt = gauss_points[GP - 1];
		real min_rho = std::numeric_limits<real>::max();
		for (integer gx = 0; gx != GP; ++gx) {
			for (integer gy = 0; gy != GP; ++gy) {
				for (integer gz = 0; gz != GP; ++gz) {
					auto u = value_at(U[rk][ii], gpt[gx], gpt[gy], gpt[gz]);
					min_rho = std::min(min_rho, u[d0i]);
				}
			}
		}
		for (integer gx = 0; gx != GP; ++gx) {
			for (integer gy = 0; gy != GP; ++gy) {
				const auto rho1 = value_at(U[rk][ii], gpt[gx], gpt[gy], +real(1))[d0i];
				const auto rho2 = value_at(U[rk][ii], gpt[gx], gpt[gy], -real(1))[d0i];
				const auto rho3 = value_at(U[rk][ii], gpt[gx], +real(1), gpt[gy])[d0i];
				const auto rho4 = value_at(U[rk][ii], gpt[gx], -real(1), gpt[gy])[d0i];
				const auto rho5 = value_at(U[rk][ii], +real(1), gpt[gx], gpt[gy])[d0i];
				const auto rho6 = value_at(U[rk][ii], -real(1), gpt[gx], gpt[gy])[d0i];
				min_rho = std::min(min_rho, std::min(std::min(rho5, rho6), std::min(std::min(rho1, rho2), std::min(rho3, rho4))));
			}
		}
		if (min_rho < rho_floor) {
			U[rk][ii][0][d0i] = std::max(rho_floor, U[rk][ii][0][d0i]);
			for (integer ppp = 1; ppp != HYDRO_PPP; ++ppp) {
				U[rk][ii][ppp] = real(0);
			}
		}
	}

}

void hydro::set_compute(bool val) {
	for (integer i = 1; i != HYDRO_NX - 1; ++i) {
		for (integer j = 1; j != HYDRO_NX - 1; ++j) {
			for (integer k = 1; k != HYDRO_NX - 1; ++k) {
				const integer ii = gindex(i, j, k);
				if (is_child_amr[ii]) {
					is_compute_cell[ii] = true;
				} else if (is_amr[ii]) {
					is_compute_cell[ii] = false;
				} else {
					is_compute_cell[ii] = val;
				}
			}
		}
	}
}

void hydro::set_amr(integer face, bool val) {
	is_amr_face[face] = val;
	integer xlb, xub, ylb, yub, zlb, zub;
	xlb = ylb = zlb = 1;
	xub = yub = zub = HYDRO_NX - 1;
	if (face == 1) {
		xlb = HYDRO_NX - 3;
	} else if (face == 0) {
		xub = 3;
	} else if (face == 3) {
		ylb = HYDRO_NX - 3;
	} else if (face == 2) {
		yub = 3;
	} else if (face == 5) {
		zlb = HYDRO_NX - 3;
	} else if (face == 4) {
		zub = 3;
	}
	for (integer i = xlb; i != xub; ++i) {
		for (integer j = ylb; j != yub; ++j) {
			for (integer k = zlb; k != zub; ++k) {
				is_amr[gindex(i, j, k)] = val;
			}
		}
	}
}

void hydro::set_child_amr(integer face, bool val) {
	is_child_amr_face[face] = val;
	integer xlb, xub, ylb, yub, zlb, zub;
	xlb = ylb = zlb = 0;
	xub = yub = zub = HYDRO_NX;
	if (face == 1) {
		xlb = HYDRO_NX - 2;
	} else if (face == 0) {
		xub = 2;
	} else if (face == 3) {
		ylb = HYDRO_NX - 2;
	} else if (face == 2) {
		yub = 2;
	} else if (face == 5) {
		zlb = HYDRO_NX - 2;
	} else if (face == 4) {
		zub = 2;
	}
	for (integer i = xlb; i != xub; ++i) {
		for (integer j = ylb; j != yub; ++j) {
			for (integer k = zlb; k != zub; ++k) {
				is_child_amr[gindex(i, j, k)] = val;
			}
		}
	}
}

gforce_t hydro::gforce_at(const std::valarray<real>& psi, const std::valarray<std::valarray<real>>& U, real x, real y, real z) const {
	std::array<real, NDIM> dist = { x, y, z };

//	dist[0] = dist[1] = dist[2] = real(0);

	gforce_t g;
	std::vector<real> this_psi(std::begin(psi), std::end(psi));
	for (integer d = 0; d != NDIM; ++d) {
		dist[d] *= dx / real(2);
	}
	exafmm.L2P(this_psi, dist, g[3], g[0], g[1], g[2]);
	gforce_t self_g = self_gforce(U, x, y, z);
	for (integer i = 0; i != NDIM; ++i) {
		g[i] += self_g[i];
	}
	return g;
}

hydro::hydro(real _dx, real _x0, real _y0, real _z0) {
	dx = _dx;
	x0 = _x0;
	y0 = _y0;
	z0 = _z0;
	psi = valarray_3d_alloc(HYDRO_RK, HYDRO_N3, FMM_PP);
	U = valarray_4d_alloc(HYDRO_RK + 1, HYDRO_N3, HYDRO_PPP, HYDRO_NF);
	dU = valarray_4d_alloc(HYDRO_RK, HYDRO_N3, HYDRO_PPP, HYDRO_NF);
	for (integer i = 0; i < HYDRO_RK + 1; ++i) {
		for (integer j = 0; j != HYDRO_N3; ++j) {
			for (integer k = 0; k != HYDRO_PPP; ++k) {
				for (integer l = 0; l != HYDRO_NF; ++l) {
					U[i][j][k][l] = real(0);
				}
			}
		}
	}
	for (integer i = 0; i < HYDRO_RK; ++i) {
		for (integer j = 0; j != HYDRO_N3; ++j) {
			for (integer k = 0; k != HYDRO_PPP; ++k) {
				for (integer l = 0; l != HYDRO_NF; ++l) {
					dU[i][j][k][l] = real(0);
				}
			}
		}
	}
	std::fill(is_amr.begin(), is_amr.end(), false);
	std::fill(is_compute_cell.begin(), is_compute_cell.end(), false);
	std::fill(is_child_amr.begin(), is_child_amr.end(), false);
	std::fill(is_amr_face.begin(), is_amr_face.end(), false);
	std::fill(is_child_amr_face.begin(), is_child_amr_face.end(), false);
}

real Ylm2Plmn[FMM_PP][HYDRO_PPP] = { { real(0) } };

struct Ylm_2_Plmn_t {
	real c;
	integer L, M;
};

struct Plmn_2_Ylm_t {
	real c;
	integer l, m, n;
};

std::list<Ylm_2_Plmn_t> Ylm_2_Plmn[HYDRO_PPP];
std::list<Plmn_2_Ylm_t> Plmn_2_Ylm[FMM_PP];

std::vector<real> hydro::get_gravity_sources(integer rk) const {
	const auto& gpt = gauss_points[HYDRO_P - 1];
	const auto& gwt = gauss_weights[HYDRO_P - 1];
	std::vector<real> src(FMM_PP * NX * NX * NX);
	const real dv = std::pow(dx, 3) / real(8);
	for (integer i = 1; i != HYDRO_NX - 1; ++i) {
		for (integer j = 1; j != HYDRO_NX - 1; ++j) {
			for (integer k = 1; k != HYDRO_NX - 1; ++k) {
				std::vector<real> this_M(FMM_PP, real(0));
				const integer iii = gindex(i, j, k);
				std::array<real, NDIM> dist;
				for (integer gx = 0; gx != HYDRO_P; ++gx) {
					const real x = gpt[gx];
					dist[0] = x * dx / real(2);
					for (integer gy = 0; gy != HYDRO_P; ++gy) {
						const real y = gpt[gy];
						dist[1] = y * dx / real(2);
						for (integer gz = 0; gz != HYDRO_P; ++gz) {
							const real z = gpt[gz];
							dist[2] = z * dx / real(2);
							const real mass = value_at(U[rk][iii], x, y, z)[d0i] * gwt[gx] * gwt[gy] * gwt[gz] * dv;
							exafmm.P2M(this_M, dist, mass);
						}
					}
				}
				for (integer L = 0; L != FMM_P; ++L) {
					for (integer M = -L; M <= L; ++M) {
						const integer PP = L * L + L + M;
						const auto yi = PP * NX * NX * NX + (i * NX + j) * NX + k - NX * NX - NX - 1;
						src[yi] = this_M[PP];
					}
				}
			}
		}
	}
	return src;
}

void hydro::set_gravity_sources(const std::vector<real>& src, integer rk) {

	for (integer i = 1; i != HYDRO_NX - 1; ++i) {
		for (integer j = 1; j != HYDRO_NX - 1; ++j) {
			for (integer k = 1; k != HYDRO_NX - 1; ++k) {
				const integer iii = gindex(i, j, k);
				for (integer L = 0; L != FMM_P; ++L) {
					for (integer M = -L; M <= L; ++M) {
						const auto ind = NX * NX * NX * (L * L + L + M) + NX * NX * i + NX * j + k - NX * NX - NX - 1;
						psi[rk][iii][L * L + L + M] = src[ind];
					}
				}
			}
		}
	}
}

void hydro_initialize() {
	const integer GP = HYDRO_P;
	const auto& gpt = gauss_points[GP - 1];
	const auto& gwt = gauss_weights[GP - 1];
	for (integer i = 0; i != FMM_PP; ++i) {
		for (integer j = 0; j != HYDRO_PPP; ++j) {
			Ylm2Plmn[i][j] = real(0);
		}
	}
	for (integer i = -1; i <= +1; ++i) {
		for (integer j = -1; j <= +1; ++j) {
			for (integer k = -1; k <= +1; ++k) {
				for (integer x1 = 0; x1 != HYDRO_P; ++x1) {
					for (integer y1 = 0; y1 != HYDRO_P; ++y1) {
						for (integer z1 = 0; z1 != HYDRO_P; ++z1) {
							for (integer x2 = 0; x2 != HYDRO_P; ++x2) {
								for (integer y2 = 0; y2 != HYDRO_P; ++y2) {
									for (integer z2 = 0; z2 != HYDRO_P; ++z2) {
										const real dx2 = std::pow((gpt[x2] - gpt[x1]) / real(2) + real(i), 2);
										const real dy2 = std::pow((gpt[y2] - gpt[y1]) / real(2) + real(j), 2);
										const real dz2 = std::pow((gpt[z2] - gpt[z1]) / real(2) + real(k), 2);
										const real r = std::sqrt(dx2 + dy2 + dz2);
										real rinv;
										if (r > real(0)) {
											rinv = real(1.0) / r;
										} else {
											rinv = real(0);
										}
										rad_inv[i + 1][j + 1][k + 1][x2][y2][z2][x1][y1][z1] = rinv;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	for (integer l = 0; l != HYDRO_P; ++l) {
		for (integer m = 0; m != HYDRO_P - l; ++m) {
			for (integer n = 0; n != HYDRO_P - l - m; ++n) {
				const integer lmn = hydro::pindex(l, m, n);
				for (integer gx = 0; gx != GP; ++gx) {
					for (integer gy = 0; gy != GP; ++gy) {
						for (integer gz = 0; gz != GP; ++gz) {
							const real x = gpt[gx];
							const real y = gpt[gy];
							const real z = gpt[gz];
							const real r = std::sqrt(x * x + y * y + z * z);
							const real theta = std::acos(z / r);
							const real phi = std::atan2(y, x);
							const real px = LegendreP(l, x) * gwt[gx];
							const real py = LegendreP(m, y) * gwt[gy];
							const real pz = LegendreP(n, z) * gwt[gz];
							std::vector<real> Ylm(FMM_PP);
							exafmm.evalMultipole(r, theta, phi, Ylm);
							for (integer L = 0; L != FMM_P; ++L) {
								if (r > real(0) || L == 0) {
									for (integer M = -L; M <= L; ++M) {
										const integer LM = L * L + L + M;
										Ylm2Plmn[LM][lmn] += Ylm[LM] * px * py * pz;
									}
								}
							}
						}

					}
				}
			}
		}
	}
	printf("!!!!!!!!\n");
	for (integer L = 0; L != FMM_P; ++L) {
		for (integer M = -L; M <= L; ++M) {
			const integer LM = L * L + L + M;
			for (integer l = 0; l != HYDRO_P; ++l) {
				for (integer m = 0; m != HYDRO_P - l; ++m) {
					for (integer n = 0; n != HYDRO_P - l - m; ++n) {
						const integer lmn = hydro::pindex(l, m, n);
						const real num = Ylm2Plmn[LM][lmn];
						if (std::abs(num) > 1.0e-10) {
							const real p = num;
							const real y = num;
							Ylm_2_Plmn_t yd;
							Plmn_2_Ylm_t pd;
							yd.c = y;
							yd.L = L;
							yd.M = M;
							pd.c = y;
							pd.l = l;
							pd.m = m;
							pd.n = n;
							Ylm_2_Plmn[lmn].push_back(yd);
							Plmn_2_Ylm[LM].push_back(pd);
							printf("%li %li %li - %li %li - %e %e\n", l, m, n, L, M, p, y);
						}
					}
				}
			}
		}
	}
}

