/*
 * exafmm.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: dmarce1
 */

#define EXAFMM_CPP
#include "exafmm.hpp"

#include <array>

#define ODDEVEN(n) real((((n) & 1) == 1) ? -1 : 1)

#define COMPLEX_MULT( ar, ai, br, bi, cr, ci) \
		(ar) = (br)*(cr)-(bi)*(ci);           \
		(ai) = (br)*(ci)+(bi)*(cr)

#define COMPLEX_MULT_ADD( ar, ai, br, bi, cr, ci) \
		(ar) += (br)*(cr)-(bi)*(ci);           \
		(ai) += (br)*(ci)+(bi)*(cr)

#define SGN(i) real((i) > 0 ? 1 : ((i)<0 ? -1 : 0))

std::array<real,FMM_P> factorial;

std::array<real,FMM_P*FMM_P> prefactor;

std::array<real,FMM_P*FMM_P> Anm;

std::array<real,FMM_P*FMM_P*FMM_P*FMM_P> Cnm_r;

std::array<real,FMM_P*FMM_P*FMM_P*FMM_P> Cnm_i;

void exafmm_kernel::M2L(std::vector<real>& CiL, const std::vector<real> CjM,
		const std::array<std::vector<real>, NDIM>& d, integer N, std::vector<real>& L_r, std::vector<real>& L_i, std::vector<real>& Ynm) {
	integer Nynm;
	Nynm = (((N - 1) / 64) + 1) * 64;
#pragma vector aligned
#pragma simd
	for (integer i = 0; i != N; ++i) {
		real rho = std::sqrt(d[0][i] * d[0][i] + d[1][i] * d[1][i] + d[2][i] * d[2][i]);
		real theta = std::acos(d[2][i] / rho);
		real phi = std::atan2(d[1][i], d[0][i]);
		real x = std::cos(theta);                              // x = cos(theta)
		real y = std::sin(theta);                              // y = sin(theta)
		real fact = 1;                                   // Initialize 2 * m + 1
		real pn = 1;                        // Initialize Legendre polynomial Pn
		real rhom = real(1.0) / rho;                          // Initialize rho^(-m-1)
#pragma novector
		for (int m = 0; m != FMM_P; ++m) {                     // Loop over m in Ynm
			real eim_r = std::cos(real(m) * phi);
			real eim_i = std::sin(real(m) * phi);
			real p = pn;                  //  Associated Legendre polynomial Pnm
			int npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
			int nmn = m * m;                          //  Index of Ynm for m < 0
			Ynm[npn * Nynm + i] = rhom * p * prefactor[npn] * eim_r; //  rho^(-m-1) * Ynm for m > 0
			if (npn != nmn) {
				Ynm[nmn * Nynm + i] = rhom * p * prefactor[npn] * eim_i; //  rho^(-m-1) * Ynm for m > 0
			}
			real p1 = p;                                              //  Pnm-1
			p = x * real(2 * m + 1) * p1;          //  Pnm using recurrence relation
			rhom /= rho;                                          //  rho^(-m-1)
			real rhon = rhom;                                     //  rho^(-n-1)
#pragma novector
			for (int n = m + 1; n != FMM_P; ++n) {            //  Loop over n in Ynm
				int npm = n * n + n + m;             //   Index of Ynm for m > 0
				int nmm = n * n + n - m;             //   Index of Ynm for m < 0
				Ynm[npm * Nynm + i] = rhon * p * prefactor[npm] * eim_r; //   rho^n * Ynm for m > 0
				if (npm != nmm) {
					Ynm[nmm * Nynm + i] = rhon * p * prefactor[npm] * eim_i; //   rho^n * Ynm for m > 0
				}
				real p2 = p1;                                         //   Pnm-2
				p1 = p;                                               //   Pnm-1
				p = (x * real(2 * n + 1) * p1 - real(n + m) * p2) / real(n - m + 1); //   Pnm using recurrence relation
				rhon /= rho;                                     //   rho^(-n-1)
			}                                         //  End loop over n in Ynm
			pn = -pn * fact * y;                                      //  Pn
			fact += real(2);                                             //  2 * m + 1
		}                                              // End loop over m in Ynm
	}

	for (integer j = 0; j != FMM_P; ++j) {
		for (integer k = 0; k <= j; ++k) {
			const integer jkp = j * j + j + k;
			const integer jkm = j * j + j - k;
#pragma vector aligned
#pragma simd
			for (integer i = 0; i != N; ++i) {
				L_r[i] = L_i[i] = real(0.0);
			}
			for (integer n = 0; n != FMM_P - j; ++n) {
				for (integer m = -n; m <= +n; ++m) {
					const integer nn = n * n + n;
					const integer nj = (n + j) * ((n + j) + 1);
					const integer jknm = jkp * FMM_P * FMM_P + n * n + n + m;
					const integer nmp = nn + std::abs(m);
					const integer nmm = nn - std::abs(m);
					const integer jnkmp = nj + std::abs(m - k);
					const integer jnkmm = nj - std::abs(m - k);
					real tmp_r, tmp_i;
					const real sgn = SGN(m-k);
					COMPLEX_MULT(tmp_r, tmp_i, CjM[nmp], SGN(m) * CjM[nmm], Cnm_r[jknm], Cnm_i[jknm]);
					const auto Yp = Ynm.data() + Nynm * jnkmp;
					const auto Ym = Ynm.data() + Nynm * jnkmm;
#pragma vector aligned
#pragma simd
					for (integer i = 0; i != N; ++i) {
						COMPLEX_MULT_ADD(L_r[i], L_i[i], tmp_r, tmp_i, Yp[i], sgn * Ym[i]);
					}
				}
			}
			auto Cp = CiL.data() + N * jkp;
			auto Cm = CiL.data() + N * jkm;
//#pragma vector aligned
#pragma simd
			for (integer i = 0; i != N; ++i) {
				Cp[i] = L_r[i];
				Cm[i] = (k == 0) ? L_r[i] : L_i[i];
			}
		}
	}

}

void exafmm_kernel::cart2sph(real& r, real& theta, real& phi, std::array<real, NDIM> dist) {
	r = std::sqrt(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);      // r = sqrt(x^2 + y^2 + z^2)
	if (r < EPS) {                                             // If r == 0
		theta = real(0);               //  theta can be anything so we set it to 0
	} else {                                                    // If r != 0
		theta = std::acos(dist[2] / r);                   //  theta = acos(z / r)
	}                                                   // End if for r == 0
	phi = std::atan2(dist[1], dist[0]);
}

void exafmm_kernel::M2M(std::vector<real>& CiM, const std::vector<real>& CjM, const std::array<real, NDIM>& dist,
		const integer N) {
	std::vector<real> Ynm(FMM_P * FMM_P);
	std::vector<real> M_r(N), M_i(N);
	real rho, theta, phi;
	cart2sph(rho, theta, phi, dist);
	evalMultipole(rho, theta, -phi, Ynm);
	for (integer j = 0; j != FMM_P; ++j) {
		for (integer k = 0; k <= j; ++k) {
			const integer jkp = j * j + j + k;
			const integer jkm = j * j + j - k;
#pragma vector aligned
#pragma simd
			for (integer i = 0; i != N; ++i) {
				M_r[i] = M_i[i] = real(0.0);
			}
			for (integer n = 0; n <= j; ++n) {
				for (integer m = std::max(n - j + k, -n); m <= std::min(j - n + k, +n); ++m) {
					const integer nn = n * n + n;
					const integer nj = (j - n) * (j - n) + j - n;
					const integer jnkm = nj + k - m;
					const integer jnpkm = nj + std::abs(k - m);
					const integer jnmkm = nj - std::abs(k - m);
					const integer nmp = nn + std::abs(m);
					const integer nmm = nn - std::abs(m);
					const auto Mj_r = CjM.data() + N * jnpkm;
					const auto Mj_i = CjM.data() + N * jnmkm;
					const real tmp = Anm[nmp] * Anm[jnkm]
							/ Anm[jkp]* ODDEVEN((std::abs(k) - std::abs(m) - std::abs(k - m)) / 2 + n);
					const real sgn_km = SGN(k-m);
					const real Y_r = tmp * Ynm[nmp];
					const real Y_i = SGN(m) * tmp * Ynm[nmm];
#pragma vector aligned
#pragma simd
					for (integer i = 0; i != N; ++i) {
						COMPLEX_MULT_ADD(M_r[i], M_i[i], Y_r, Y_i, Mj_r[i], sgn_km * Mj_i[i]);
					}
				}
			}
			auto Mi_r = CiM.data() + N * jkp;
			auto Mi_i = CiM.data() + N * jkm;
#pragma vector aligned
#pragma simd
			for (integer i = 0; i != N; ++i) {
				Mi_r[i] += M_r[i];
				Mi_i[i] += (jkm == jkp) ? 0.0 : M_i[i];
			}
		}
	}
}

void exafmm_kernel::L2L(std::vector<real>& CiL, const std::vector<real>& CjL, const std::array<real, NDIM>& dist,
		const integer N) {
	std::vector<real> Ynm(FMM_P * FMM_P);
	real rho, theta, phi;
	std::vector<real> L_r(N), L_i(N);
	cart2sph(rho, theta, phi, dist);
	evalMultipole(rho, theta, phi, Ynm);
	for (integer j = 0; j != FMM_P; ++j) {
		for (integer k = 0; k <= j; ++k) {
			integer jkp = j * j + j + k;
			integer jkm = j * j + j - k;
#pragma vector aligned
#pragma simd
			for (integer i = 0; i != N; ++i) {
				L_r[i] = L_i[i] = 0.0;
			}
			for (integer n = j; n != FMM_P; ++n) {
				for (integer m = j - n + k; m <= n - j + k; ++m) {
					const integer nn = n * n + n;
					const integer nj = (n - j) * ((n - j) + 1);
					const integer npm = nn + std::abs(m);
					const integer nmm = nn - std::abs(m);
					const integer jnpkm = nj + std::abs(m - k);
					const integer jnmkm = nj - std::abs(m - k);
					const auto Lj_r = CjL.data() + N * npm;
					const auto Lj_i = CjL.data() + N * nmm;
					const real sgn = SGN(m);
					real tmp = std::pow(-real(1.0), real(std::abs(m) - std::abs(k) - std::abs(m - k)) / 2) * Anm[jnpkm] * Anm[jkp]
							/ Anm[npm];
					const real Y_r = Ynm[jnpkm] * tmp;
					const real Y_i = SGN(m-k) * Ynm[jnmkm] * tmp;
#pragma vector aligned
#pragma simd
					for (integer i = 0; i != N; ++i) {
						COMPLEX_MULT_ADD(L_r[i], L_i[i], Y_r, Y_i, Lj_r[i], sgn * Lj_i[i]);
					}
				}
			}
			auto Li_r = CiL.data() + N * jkp;
			auto Li_i = CiL.data() + N * jkm;
#pragma vector aligned
#pragma simd
			for (integer i = 0; i != N; ++i) {
				Li_r[i] = L_r[i];
				Li_i[i] = (k == 0) ? L_r[i] : L_i[i];
			}
		}
	}
}

void exafmm_kernel::evalMultipole(real rho, real theta, real phi, std::vector<real>& Ynm) {
	real x = std::cos(theta);                              // x = cos(theta)
	real y = std::sin(theta);                              // y = sin(theta)
	real fact = 1;                                   // Initialize 2 * m + 1
	real pn = 1;                        // Initialize Legendre polynomial Pn
	real rhom = 1;                                       // Initialize rho^m
	for (int m = 0; m != FMM_P; ++m) {                     // Loop over m in Ynm
		real eim_r = std::cos(real(m) * phi);
		real eim_i = std::sin(real(m) * phi);
		real p = pn;                  //  Associated Legendre polynomial Pnm
		int npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
		int nmn = m * m;                          //  Index of Ynm for m < 0
		Ynm[npn] = rhom * p * prefactor[npn] * eim_r; //  rho^m * Ynm for m > 0
		if (npn != nmn) {
			Ynm[nmn] = rhom * p * prefactor[npn] * eim_i; //  rho^m * Ynm for m > 0
		}
		real p1 = p;                                              //  Pnm-1
		p = x *real(2 * m + 1) * p1;          //  Pnm using recurrence relation
		rhom *= rho;                                              //  rho^m
		real rhon = rhom;                                         //  rho^n
		for (int n = m + 1; n != FMM_P; ++n) {            //  Loop over n in Ynm
			int npm = n * n + n + m;             //   Index of Ynm for m > 0
			int nmm = n * n + n - m;             //   Index of Ynm for m < 0
			Ynm[npm] = rhon * p * prefactor[npm] * eim_r;     //   rho^n * Ynm
			if (npm != nmm) {
				Ynm[nmm] = rhon * p * prefactor[npm] * eim_i;     //   rho^n * Ynm
			}
			real p2 = p1;                                         //   Pnm-2
			p1 = p;                                               //   Pnm-1
			p = (x * real(2 * n + 1) * p1 - (n + m) * p2) / real(n - m + 1); //   Pnm using recurrence relation
			rhon *= rho;                                   //   Update rho^n
		}                                         //  End loop over n in Ynm
		pn = -pn * fact * y;                                      //  Pn
		fact += real(2);                                             //  2 * m + 1
	}                                              // End loop over m in Ynm
}

exafmm_kernel::exafmm_kernel() {
	const complex I(0., 1.);                               // Imaginary unit
	if( FMM_P == 0 ) {
		return;
	}
	factorial[0] = 1;                                // Initialize factorial
	for (int n = 1; n < FMM_P; ++n) {                              // Loop to P
		factorial[n] = factorial[n - 1] * n;                        //  n!
	}                                                       // End loop to P

	for (int n = 0; n != FMM_P; ++n) {                     // Loop over n in Anm
		for (int m = -n; m <= n; ++m) {               //  Loop over m in Anm
			int nm = n * n + n + m;                        //   Index of Anm
			int nabsm = abs(m);                                     //   |m|
			real fnmm = real(1.0);                        //   Initialize (n - m)!
			for (int i = 1; i <= n - m; ++i)
				fnmm *= i;                  //   (n - m)!
			real fnpm = real(1.0);                        //   Initialize (n + m)!
			for (int i = 1; i <= n + m; ++i)
				fnpm *= i;                  //   (n + m)!
			real fnma = real(1.0);                      //   Initialize (n - |m|)!
			for (int i = 1; i <= n - nabsm; ++i)
				fnma *= i;              //   (n - |m|)!
			real fnpa = real(1.0);                      //   Initialize (n + |m|)!
			for (int i = 1; i <= n + nabsm; ++i)
				fnpa *= i;              //   (n + |m|)!
			prefactor[nm] = std::sqrt(fnma / fnpa); //   sqrt( (n - |m|)! / (n + |m|)! )
			Anm[nm] = ODDEVEN(n) / std::sqrt(fnmm * fnpm); //   (-1)^n / sqrt( (n + m)! / (n - m)! )
		}                                         //  End loop over m in Anm
	}                                              // End loop over n in Anm

	for (int j = 0, jk = 0, jknm = 0; j != FMM_P; ++j) { // Loop over j in Cjknm
		for (int k = -j; k <= j; ++k, ++jk) {       //  Loop over k in Cjknm
			for (int n = 0, nm = 0; n != FMM_P; ++n) { //   Loop over n in Cjknm
				for (int m = -n; m <= n; ++m, ++nm, ++jknm) { //    Loop over m in Cjknm
					if (j + n < FMM_P) {
						const int jnkm = (j + n) * (j + n) + j + n + m - k; //     Index C_{j+n}^{m-k}
						auto tmp = std::pow(I, real(abs(k - m) - abs(k) - abs(m))) //     Cjknm
						* real(ODDEVEN(j) * Anm[nm] * Anm[jk] / Anm[jnkm]);
						Cnm_r[jknm] = tmp.real();
						Cnm_i[jknm] = tmp.imag();

					}
				}
			}
		}
	}
}
