/*
 * hydro.cpp
 *
 *  Created on: Mar 6, 2015
 *      Author: dmarce1
 */

#include "hydro.hpp"
#include "math.hpp"
#include "rk.hpp"

const integer di = HYDRO_NX * HYDRO_NX;
const integer dj = HYDRO_NX;
const integer dk = 1;

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

std::function<state(real, real, real)> init_func = sod_shock;

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

	for (integer i = 0; i != HYDRO_NX; ++i) {
		for (integer j = 0; j != HYDRO_NX; ++j) {
			for (integer k = 0; k != HYDRO_NX; ++k) {
				integer ii = gindex(i, j, k);
				real x, y, z;
				x = x0 + (i - 0.5) * dx;
				y = y0 + (j - 0.5) * dx;
				z = z0 + (k - 0.5) * dx;
				for (integer pp = 0; pp != HYDRO_PPP; ++pp) {
					U[0][ii][pp] = real(0);
				}
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

/*** Change sample points***********/
std::vector<real> hydro::output_data() const {
	const auto& gpt = gauss_points[HYDRO_P - 1];
	std::vector<real> data(SILO_N3 * HYDRO_NF);
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
				for (integer f = 0; f != HYDRO_NF; ++f) {
					*iter++ = u[f];
				}
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
	const integer NGF = std::max(HYDRO_P - 1, integer(1));
	const integer NGC = std::max(HYDRO_P - 1, integer(1));
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
						const state& urx = value_at(U[rk][i], -real(1), x1, x2);
						const state& ury = value_at(U[rk][i], x1, -real(1), x2);
						const state& urz = value_at(U[rk][i], x1, x2, -real(1));
						const state& ulx = value_at(U[rk][i - di], +real(1), x1, x2);
						const state& uly = value_at(U[rk][i - dj], x1, +real(1), x2);
						const state& ulz = value_at(U[rk][i - dk], x1, x2, +real(1));
						state fx = Riemann_flux(urx, ulx, 0);
						state fy = Riemann_flux(ury, uly, 1);
						state fz = Riemann_flux(urz, ulz, 2);
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
								const state& u = value_at(U[rk][i], x, y, z);
								state fx = flux(u, 0);
								state fy = flux(u, 1);
								state fz = flux(u, 2);
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
	return dx / amax;
}

std::vector<real> hydro::flux_correct_pack(integer rk) const {
	const integer NGF = std::max(HYDRO_P - 1, integer(1));
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
									const state& urx = value_at(U[rk][i], -real(1), gfpt[g1], gfpt[g2]);
									const state& ulx = value_at(U[rk][i - di], +real(1), gfpt[g1], gfpt[g2]);
									state fx = Riemann_flux(urx, ulx, 0);
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
									const state& ury = value_at(U[rk][i], gfpt[g1], -real(1), gfpt[g2]);
									const state& uly = value_at(U[rk][i - dj], gfpt[g1], +real(1), gfpt[g2]);
									state fy = Riemann_flux(ury, uly, 1);
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
									const state& urz = value_at(U[rk][i], gfpt[g1], gfpt[g2], -real(1));
									const state& ulz = value_at(U[rk][i - dk], gfpt[g1], gfpt[g2], +real(1));
									state fz = Riemann_flux(urz, ulz, 2);
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
	for (integer l = 0; l != HYDRO_P; ++l) {
		for (integer m = 0; m != HYDRO_P - l; ++m) {
			for (integer n = 0; n != HYDRO_P - l - m; ++n) {
				for (integer j = 0; j != HYDRO_NX; ++j) {
					for (integer k = 0; k != HYDRO_NX; ++k) {
						const integer p = pindex(l, m, n);
						if (dir == 0) {
							U[rk][gindex(0, j, k)][p] = U[rk][gindex(1, j, k)][p];
						} else if (dir == 1) {
							U[rk][gindex(HYDRO_NX - 1, j, k)][p] = U[rk][gindex(HYDRO_NX - 2, j, k)][p];
						} else if (dir == 2) {
							U[rk][gindex(j, 0, k)][p] = U[rk][gindex(j, 1, k)][p];
						} else if (dir == 3) {
							U[rk][gindex(j, HYDRO_NX - 1, k)][p] = U[rk][gindex(j, HYDRO_NX - 2, k)][p];
						} else if (dir == 4) {
							U[rk][gindex(j, k, 0)][p] = U[rk][gindex(j, k, 1)][p];
						} else if (dir == 5) {
							U[rk][gindex(j, k, HYDRO_NX - 1)][p] = U[rk][gindex(j, k, HYDRO_NX - 2)][p];
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

std::valarray<real> minmod(const std::valarray<real>& a, const std::valarray<real>& b) {
	state c(HYDRO_NF);
	constexpr auto half = real(1) / real(2);
	for (integer i = 0; i != HYDRO_NF; ++i) {
		c[i] = (std::copysign(half, a[i]) + std::copysign(half, b[i])) * std::min(std::abs(a[i]), std::abs(b[i]));
	}
	return c;
}
;

bool real_eq(real a, real b) {
	constexpr real delta = 1.0e-11;
	return std::abs(a - b) / (std::abs(a) + std::abs(b)) < delta;
}

void hydro::apply_limiter(integer rk) {

	auto& u = U[rk];
	auto ux = u;
	auto uy = u;
	auto uz = u;
	auto utmp = u;

	for (integer i = di; i < HYDRO_N3 - di; ++i) {
		for (integer m = 0; m < HYDRO_P - 1; ++m) {
			for (integer n = 0; n < HYDRO_P - m - 1; ++n) {
				for (integer l = HYDRO_P - 2 - n - m; l >= 0; --l) {
					const integer p = pindex(l, m, n);
					const integer pp1 = pindex(l + 1, m, n);
					state dup = (u[i + di][p] - u[i][p]) / real(2 * l + 1);
					state dum = (u[i][p] - u[i - di][p]) / real(2 * l + 1);
					state du0 = u[i][pp1];
					dup = con_to_characteristic(u[i][0], dup, 0);
					dum = con_to_characteristic(u[i][0], dum, 0);
					du0 = con_to_characteristic(u[i][0], du0, 0);
					utmp[i][pp1] = characteristic_to_con(u[i][0], minmod(du0, minmod(dup, dum)), 0);
				}
				for (integer f = 0; f != HYDRO_NF; ++f) {
					for (integer l = HYDRO_P - 2 - n - m; l >= 0; --l) {
						const integer pp1 = pindex(l + 1, m, n);
						if (!real_eq(utmp[i][pp1][f], u[i][pp1][f])) {
							ux[i][pp1][f] = utmp[i][pp1][f];
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
					dup = con_to_characteristic(u[i][0], dup, 1);
					dum = con_to_characteristic(u[i][0], dum, 1);
					du0 = con_to_characteristic(u[i][0], du0, 1);
					utmp[i][pp1] = characteristic_to_con(u[i][0], minmod(du0, minmod(dup, dum)), 1);
				}
				for (integer f = 0; f != HYDRO_NF; ++f) {
					for (integer m = HYDRO_P - 2 - l - n; m >= 0; --m) {
						const integer pp1 = pindex(l, m + 1, n);
						if (!real_eq(utmp[i][pp1][f], u[i][pp1][f])) {
							uy[i][pp1][f] = utmp[i][pp1][f];
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
					dup = con_to_characteristic(u[i][0], dup, 2);
					dum = con_to_characteristic(u[i][0], dum, 2);
					du0 = con_to_characteristic(u[i][0], du0, 2);
					utmp[i][pp1] = characteristic_to_con(u[i][0], minmod(du0, minmod(dup, dum)), 2);
				}
				for (integer f = 0; f != HYDRO_NF; ++f) {
					for (integer n = HYDRO_P - 2 - m - l; n >= 0; --n) {
						const integer pp1 = pindex(l, m, n + 1);
						if (!real_eq(utmp[i][pp1][f], u[i][pp1][f])) {
							uz[i][pp1][f] = utmp[i][pp1][f];
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
			u[i][pp] = minmod(uz[i][pp], minmod(ux[i][pp], uy[i][pp]));
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

hydro::hydro(real _dx, real _x0, real _y0, real _z0) {
	dx = _dx;
	x0 = _x0;
	y0 = _y0;
	z0 = _z0;
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
	initialize();
}
