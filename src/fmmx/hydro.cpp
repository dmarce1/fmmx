/*
 * hydro.cpp
 *
 *  Created on: Feb 4, 2015
 *      Author: dmarce1
 */

#include "hydro.hpp"
#include "lane_emden.hpp"

constexpr real hydro_vars::ro_floor;

real minmod(real a, real b) {
	if (a * b >= real(0)) {
		return std::copysign(std::min(std::abs(a), std::abs(b)), a);
	} else {
		return real(0);
	}
}

integer hydro_vars::ind3d(integer j, integer k, integer l, integer stride) {
	return l + stride * (k + stride * j);
}

void hydro_vars::store() {
	U0 = U;
	for (integer d = 0; d != NDIM; ++d) {
		Ux0[d] = Ux[d];
	}
}

real hydro_vars::cell_mass(integer i) const {
	return U[d0_i][i] * dx * dx * dx;
}

void hydro_vars::prolong_pack(integer face, integer child, std::vector<real>::iterator i) const {
	integer lb[NDIM], ub[NDIM];
	lb[0] = ((child >> 0) & 1) * NX / 2;
	lb[1] = ((child >> 1) & 1) * NX / 2;
	lb[2] = ((child >> 2) & 1) * NX / 2;
	ub[0] = lb[0] + NX / 2;
	ub[1] = lb[1] + NX / 2;
	ub[2] = lb[2] + NX / 2;
	lb[face / 2] += (face % 2) * (NX / 2 - 1);
	ub[face / 2] = lb[face / 2] + 1;
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer j = lb[0]; j != ub[0]; ++j) {
			for (integer k = lb[1]; k != ub[1]; ++k) {
				for (integer l = lb[2]; l != ub[2]; ++l) {
					*i++ = U[f][ind3d(j, k, l)];
					for (integer d = 0; d != NDIM; ++d) {
						*i++ = Ux[d][f][ind3d(j, k, l)];
					}
				}
			}
		}
	}
}

void hydro_vars::prolong_unpack(integer face, std::vector<real>::iterator i) {
	integer lb[NDIM] = { 0, 0, 0 };
	integer ub[NDIM] = { NX / 2, NX / 2, NX / 2 };
	lb[face / 2] = (face % 2) * (NX / 2 - 1);
	ub[face / 2] = lb[face / 2] + 1;
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer j = lb[0]; j != ub[0]; ++j) {
			for (integer k = lb[1]; k != ub[1]; ++k) {
				for (integer l = lb[2]; l != ub[2]; ++l) {
					real v = *i++;
					real slp[NDIM];
					for (integer j0 = 0; j0 != 2; ++j0) {
						for (integer k0 = 0; k0 != 2; ++k0) {
							for (integer l0 = 0; l0 != 2; ++l0) {
								U[f][ind3d(2 * j + j0, 2 * k + k0, 2 * l + l0)] = v;
							}
						}
					}
					for (integer d = 0; d != NDIM; ++d) {
						real v = *i++;
						slp[d] = v;
						for (integer j0 = 0; j0 != 2; ++j0) {
							for (integer k0 = 0; k0 != 2; ++k0) {
								for (integer l0 = 0; l0 != 2; ++l0) {
									Ux[d][f][ind3d(2 * j + j0, 2 * k + k0, 2 * l + l0)] = v;
								}
							}
						}
					}
				}
			}
		}
	}

}

void hydro_vars::restrict_pack(std::vector<real>& data, std::vector<bool> amr_dirs) const {
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer j = 0; j != NX; ++j) {
			for (integer k = 0; k != NX; ++k) {
				for (integer l = 0; l != NX; ++l) {
					data[f * (N3 / 8) + ind3d(j / 2, k / 2, l / 2, NX / 2)] += dU[f][ind3d(j, k, l)] * real(1.0 / 8.0);
					for (integer d = 0; d != NDIM; ++d) {
						integer offset = (d + 1) * NX * NX * NX / 8 * nf_hydro;
						data[offset + f * (N3 / 8) + ind3d(j / 2, k / 2, l / 2, NX / 2)] += dUx[d][f][ind3d(j, k, l)]
								* real(1.0 / 8.0);
					}
				}
			}
		}
	}
	integer lb[NDIM], ub[NDIM], off[NDIM];
	for (integer d = 0; d != 2 * NDIM; ++d) {
		if (amr_dirs[d]) {
			off[0] = off[1] = off[2] = 0;
			lb[0] = lb[1] = lb[2] = 0;
			ub[0] = ub[1] = ub[2] = NX;
			if (d % 2 == 0) {
				lb[d / 2] = 2;
				off[d / 2] = 2;
			} else {
				lb[d / 2] = NX - 2;
			}
			ub[d / 2] = lb[d / 2] + 1;
			for (integer f = 0; f != nf_hydro; ++f) {
				for (integer j = lb[0]; j < ub[0]; j += 2) {
					for (integer k = lb[1]; k < ub[1]; k += 2) {
						for (integer l = lb[2]; l < ub[2]; l += 2) {
							data[f * (N3 / 8) + ind3d((j - off[0]) / 2, (k - off[1]) / 2, (l - off[2]) / 2, NX / 2)] =
									real(0);
						}
					}
				}
				for (integer j = lb[0]; j != ub[0]; ++j) {
					for (integer k = lb[1]; k != ub[1]; ++k) {
						for (integer l = lb[2]; l != ub[2]; ++l) {
							data[f * (N3 / 8) + ind3d((j - off[0]) / 2, (k - off[1]) / 2, (l - off[2]) / 2, NX / 2)] +=
									FaceFlux[d / 2][f][ind3d(j, k, l, NX + 1)] * real(1.0 / 4.0);
						}
					}
				}
			}
		}
	}
}

void hydro_vars::restrict_unpack(std::vector<real>& data, integer child, std::vector<bool> amr_directions) {

	integer lb[NDIM], ub[NDIM], off[NDIM], off2[NDIM];

	for (integer d = 0; d != 2 * NDIM; ++d) {
		if (amr_directions[d]) {
			real sign;
			lb[0] = ((child >> 0) & 1) * NX / 2;
			lb[1] = ((child >> 1) & 1) * NX / 2;
			lb[2] = ((child >> 2) & 1) * NX / 2;
			off2[0] = off2[1] = off2[2] = 0;
			off[0] = lb[0];
			off[1] = lb[1];
			off[2] = lb[2];
			ub[0] = lb[0] + NX / 2;
			ub[1] = lb[1] + NX / 2;
			ub[2] = lb[2] + NX / 2;
			for (integer e = 0; e != 2 * NDIM; e += 2) {
				if (amr_directions[e]) {
					lb[e / 2]++;
				}
				if (amr_directions[e + 1]) {
					ub[e / 2]--;
				}
			}
			if (d % 2 == 0) {
				sign = real(1);
				lb[d / 2] = 1;
				off[d / 2] = 1;
				off2[d / 2]--;
			} else {
				off[d / 2] = NX / 2;
				sign = real(-1);
				lb[d / 2] = NX - 1;
			}
			ub[d / 2] = lb[d / 2] + 1;
			for (integer f = 0; f != nf_hydro; ++f) {
				for (integer j = lb[0]; j != ub[0]; ++j) {
					for (integer k = lb[1]; k != ub[1]; ++k) {
						for (integer l = lb[2]; l != ub[2]; ++l) {
							const auto fine_flux =
									data[f * (N3 / 8) + ind3d(j - off[0], k - off[1], l - off[2], NX / 2)];
							const auto crse_flux = FaceFlux[d / 2][f][ind3d(j, k, l, NX + 1)];
							const auto dif = fine_flux - crse_flux;
							//					if (fine_flux != 0.0 || crse_flux != 0.0)
							//						printf("%e %e\n", fine_flux, crse_flux);
							dU[f][ind3d(j + off2[0], k + off2[1], l + off2[2])] -= sign * dif / dx;
						}
					}
				}
			}
		}
	}

	off[0] = lb[0] = ((child >> 0) & 1) * NX / 2;
	off[1] = lb[1] = ((child >> 1) & 1) * NX / 2;
	off[2] = lb[2] = ((child >> 2) & 1) * NX / 2;
	ub[0] = lb[0] + NX / 2;
	ub[1] = lb[1] + NX / 2;
	ub[2] = lb[2] + NX / 2;
	for (integer i = 0; i != 2 * NDIM; i += 2) {
		if (amr_directions[i]) {
			lb[i / 2]++;
		}
		if (amr_directions[i + 1]) {
			ub[i / 2]--;
		}
	}
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer j = lb[0]; j != ub[0]; ++j) {
			for (integer k = lb[1]; k != ub[1]; ++k) {
				for (integer l = lb[2]; l != ub[2]; ++l) {
					dU[f][ind3d(j, k, l)] = data[f * (N3 / 8) + ind3d(j - off[0], k - off[1], l - off[2], NX / 2)];
					for (integer d = 0; d != NDIM; ++d) {
						integer offset = (d + 1) * NX * NX * NX / 8 * nf_hydro;
						dUx[d][f][ind3d(j, k, l)] = data[offset + f * (N3 / 8)
								+ ind3d(j - off[0], k - off[1], l - off[2], NX / 2)];
					}
				}
			}
		}
	}
}

real hydro_vars::compute_du() {
	const auto dxinv = real(1) / dx;
	std::vector<std::vector<real>> F(nf_hydro, std::vector<real>(N3F, real(0)));
	for (integer f = 0; f != nf_hydro; ++f) {
		for (integer i = 0; i != N3; ++i) {
			dU[f][i] = real(0);
			for (integer d = 0; d != NDIM; ++d) {
				dUx[d][f][i] = real(0);
			}
		}
	}
	std::vector<std::vector<real>> flx(nf_hydro, std::vector<real>(N3F, real(0)));
	for (integer d = 0; d != NDIM; ++d) {
		const integer dn[3] = { d == 0 ? 1 : 0, d == 1 ? 1 : 0, d == 2 ? 1 : 0 };
		for (integer f = 0; f != nf_hydro; ++f) {
			for (integer j = 0; j != NX; ++j) {
				for (integer k = 0; k != NX; ++k) {
					for (integer l = 0; l != NX; ++l) {
						const auto iu0 = ind3d(j, k, l);
						const auto if0 = ind3d(j, k, l, NX + 1);
						const auto ifp = ind3d(j + dn[0], k + dn[1], l + dn[2], NX + 1);
						const auto Fp = FaceFlux[d][f][ifp];
						const auto Fm = FaceFlux[d][f][if0];
						const auto Fc = CenterFlux[d][f][iu0];
						const auto Xp = (Xf[d][ifp] - Xc[d][iu0]) * dxinv;
						const auto Xm = (Xf[d][if0] - Xc[d][iu0]) * dxinv;
						dU[f][iu0] += -(Fp - Fm) * dxinv;
						dUx[d][f][iu0] += real(12) * (-(Xp * Fp - Xm * Fm) + Fc) * dxinv;
					}
				}
			}
		}
	}
	return cfl_factor * dx / amax;
}

hydro_vars::hydro_vars() {
	amax = real(0);
	U = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
	dU = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
	U0 = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
	for (integer d = 0; d != NDIM; ++d) {
		Ux[d] = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
		dUx[d] = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
		Ux0[d] = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
		Xc[d] = std::vector<real>(N3);
		Xf[d] = std::vector<real>(N3F);
	}
	r = std::vector<real>(N3);
	for (integer d = 0; d != NDIM; ++d) {
		FaceFlux[d] = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3F, real(0.0)));
		CenterFlux[d] = std::vector<std::vector<real>>(nf_hydro, std::vector<real>(N3, real(0.0)));
	}

}

void hydro_vars::dump_data(std::vector<double>::iterator& l, integer index) const {
	for (integer f = 0; f != nf_hydro; ++f) {
		*l++ = double(U[f][index]);
	}
	for (integer d = 0; d != NDIM; ++d) {
		for (integer f = 0; f != nf_hydro; ++f) {
			*l++ = double(Ux[d][f][index]);
		}
	}
}

void hydro_vars::initialize(real xcorner, real ycorner, real zcorner, real _dx) {
	dx = _dx;
	for (integer i = 0; i != N3; ++i) {
		auto const j = i / (NX * NX);
		auto const k = (i / NX) % NX;
		auto const l = i % NX;
		Xc[0][i] = xcorner + (real(j) + real(0.5)) * dx - real(0.5);
		Xc[1][i] = ycorner + (real(k) + real(0.5)) * dx - real(0.5);
		Xc[2][i] = zcorner + (real(l) + real(0.5)) * dx - real(0.5);
		r[i] = std::sqrt(Xc[0][i] * Xc[0][i] + Xc[1][i] * Xc[1][i] + Xc[2][i] * Xc[2][i]);
		for (integer j = 0; j != NF; ++j) {
			Ux[0][i][j] = Ux[1][i][j] = Ux[2][i][j] = 0.0;
		}
	}
	for (integer i = 0; i != N3F; ++i) {
		auto const j = i / ((NX + 1) * (NX + 1));
		auto const k = (i / (NX + 1)) % (NX + 1);
		auto const l = i % (NX + 1);
		Xf[0][i] = xcorner + real(j) * dx - real(0.5);
		Xf[1][i] = ycorner + real(k) * dx - real(0.5);
		Xf[2][i] = zcorner + real(l) * dx - real(0.5);
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

		if (Xc[2][i] < 0.) {
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
	for (integer i = 0; i != N3; ++i) {
		if (std::abs(r[i]) < 0.125) {
			return true;
		}
	}
	return false;
}

void hydro_vars::update(real dt, integer rk, const std::vector<real>& phi, const std::vector<real>& gx,
		const std::vector<real>& gy, const std::vector<real>& gz) {
	real beta[2] = { 1.0, 0.5 };
	amax = real(0);
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
			for (integer d = 0; d != NDIM; ++d) {
				real u1 = Ux[d][f][i] + dUx[d][f][i] * dt;
				real u0 = Ux0[d][f][i];
				Ux[d][f][i] = u1 * beta[rk] + u0 * (1.0 - beta[rk]);
			}
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
		const real ek = real(0.5) * (sx * sx + sy * sy + sz * sz) / ro;
		V[pr_i][i] = std::max((fgamma - 1.0) * (et - ek), 0.0);
	}
	std::vector<real> vr(nf_hydro), vl(nf_hydro);
	for (integer j = 2 * dn[0]; j != xbdim[d] - dn[0]; ++j) {
		for (integer k = 2 * dn[1]; k != ybdim[d] - dn[1]; ++k) {
			for (integer l = 2 * dn[2]; l != zbdim[d] - dn[2]; ++l) {
				for (integer f = 0; f != nf_hydro; ++f) {
					const auto ri = ind3d(j + xoff[d], k + yoff[d], l + zoff[d], NX + 1);
					const auto vi0 = this_ind3d(j, k, l);
					const auto vim = this_ind3d(j - dn[0], k - dn[1], l - dn[2]);
					vr[f] = V[f][vi0];
					vl[f] = V[f][vim];
				}
				const auto fi = ind3d(j + xoff[d], k + yoff[d], l + zoff[d], NX + 1);
				const auto ii = ind3d(j + xoff[d], k + yoff[d], l + zoff[d]);
				const real rho_r = vr[ro_i];
				const real rho_l = vl[ro_i];
				const real v_r = vr[v0_i + d / 2];
				const real v_l = vl[v0_i + d / 2];
				const real vx_r = vr[vx_i];
				const real vx_l = vl[vx_i];
				const real vy_r = vr[vy_i];
				const real vy_l = vl[vy_i];
				const real vz_r = vr[vz_i];
				const real vz_l = vl[vz_i];
				const real pr_r = vr[pr_i];
				const real pr_l = vl[pr_i];
				const real et_r = pr_r / (fgamma - 1.0) + 0.5 * (vx_r * vx_r + vy_r * vy_r + vz_r * vz_r) * rho_r;
				const real et_l = pr_l / (fgamma - 1.0) + 0.5 * (vx_l * vx_l + vy_l * vy_l + vz_l * vz_l) * rho_l;
				const real a_r = std::sqrt(fgamma * std::max(pr_r, real(0)) / rho_r) + std::abs(v_r);
				const real a_l = std::sqrt(fgamma * std::max(pr_l, real(0)) / rho_l) + std::abs(v_l);
				const real a = std::max(a_r, a_l);
				amax = std::max(amax, a);
				FaceFlux[d / 2][d0_i][fi] = real(0.5) * (rho_r * v_r + rho_l * v_l);
				FaceFlux[d / 2][sx_i][fi] = real(0.5) * (rho_r * vx_r * v_r + rho_l * vx_l * v_l);
				FaceFlux[d / 2][sy_i][fi] = real(0.5) * (rho_r * vy_r * v_r + rho_l * vy_l * v_l);
				FaceFlux[d / 2][sz_i][fi] = real(0.5) * (rho_r * vz_r * v_r + rho_l * vz_l * v_l);
				FaceFlux[d / 2][s0_i + d / 2][fi] += real(0.5) * (pr_r + pr_l);
				FaceFlux[d / 2][et_i][fi] = real(0.5) * ((et_r + pr_r) * v_r + (et_l + pr_l) * v_l);
				FaceFlux[d / 2][d0_i][fi] -= real(0.5) * (rho_r - rho_l) * a;
				FaceFlux[d / 2][sx_i][fi] -= real(0.5) * (vx_r * rho_r - vx_l * rho_l) * a;
				FaceFlux[d / 2][sy_i][fi] -= real(0.5) * (vy_r * rho_r - vy_l * rho_l) * a;
				FaceFlux[d / 2][sz_i][fi] -= real(0.5) * (vz_r * rho_r - vz_l * rho_l) * a;
				FaceFlux[d / 2][et_i][fi] -= real(0.5) * (et_r - et_l) * a;
			}
		}
	}
	std::vector<real> vc(nf_hydro);
	for (integer j = 2 * dn[0]; j != xbdim[d] - 2 * dn[0]; ++j) {
		for (integer k = 2 * dn[1]; k != ybdim[d] - 2 * dn[1]; ++k) {
			for (integer l = 2 * dn[2]; l != zbdim[d] - 2 * dn[2]; ++l) {
				for (integer f = 0; f != nf_hydro; ++f) {
					const auto vi0 = this_ind3d(j, k, l);
					vc[f] = V[f][vi0];
				}
				const auto fi = ind3d(j + xoff[d], k + yoff[d], l + zoff[d]);
				const real rho_c = vc[ro_i];
				const real v_c = vc[v0_i + d / 2];
				const real vx_c = vc[vx_i];
				const real vy_c = vc[vy_i];
				const real vz_c = vc[vz_i];
				const real pr_c = vc[pr_i];
				const real et_c = pr_c / (fgamma - 1.0) + 0.5 * (vx_c * vx_c + vy_c * vy_c + vz_c * vz_c) * rho_c;
				CenterFlux[d / 2][d0_i][fi] = (rho_c * v_c);
				CenterFlux[d / 2][sx_i][fi] = (rho_c * vx_c * v_c);
				CenterFlux[d / 2][sy_i][fi] = (rho_c * vy_c * v_c);
				CenterFlux[d / 2][sz_i][fi] = (rho_c * vz_c * v_c);
				CenterFlux[d / 2][s0_i + d / 2][fi] += (pr_c);
				CenterFlux[d / 2][et_i][fi] = ((et_c + pr_c) * v_c);
			}
		}
	}
}
