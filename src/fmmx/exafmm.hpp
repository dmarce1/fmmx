/****
 * THIS FILE IS ADAPTED FROM EXAFMM (SEE COPYRIGHT NOTICE BELOW). THE INTERFACE HAS BEEN ALTERED
 * SLIGHTLY
 */

/*
 Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#include "defs.hpp"

#ifndef kernel_h
#define kernel_h

#include <complex>
#include <vector>
#include <array>

#define EPS2 0.0

using complex = std::complex<real>;

#define EPS (1.0e-12)

class exafmm_kernel {
public:
	static void P2M(std::vector<real>& M, const std::array<real,NDIM> dist, real m);
	static void L2P(std::vector<real>& L, const std::array<real,NDIM> dist, real&, real&, real&, real&);

	static void M2L(std::vector<real>& CiL, const std::vector<real> CjM,
			const std::array<std::vector<real>,NDIM>& dist, integer N, std::vector<real>& L_r, std::vector<real>& L_i, std::vector<real>& Ynm);

	static void cart2sph(real& r, real& theta, real& phi, std::array<real, NDIM> dist);
	static void M2M(std::vector<real>& CiM, const std::vector<real>& CjM, const std::array<real, NDIM>& dist,
			const integer N);
	static void L2L(std::vector<real>& CiL, const std::vector<real>& CjL, const std::array<real, NDIM>& dist,
			const integer N);
	static void evalMultipole(real rho, real theta, real phi, std::vector<real>& Ynm);
	static void evalMultipole(real rho, real theta, real phi, std::vector<real>& Ynm, std::vector<real>& Ynm_theta);

public:
	exafmm_kernel();
};

#endif
