/*
 * lane_emden.cpp
 *
 *  Created on: Feb 6, 2015
 *      Author: dmarce1
 */


#include "lane_emden.hpp"

static constexpr real n = 1.5;

static real fy(real y, real z, real r) {
	return z;
}

static real fz(real y, real z, real r) {
	if (r != 0.0) {
		return -(std::pow(y, n) + 2.0 * z / r);
	} else {
		return -3.0;
	}
}

real lane_emden(real r0, double dr) {
	int N;
	real dy1, dz1, y, z, r, dy2, dz2, dy3, dz3, dy4, dz4, y0, z0;
	int done = 0;
	y = 1.0;
	z = 0.0;
	N = (int) (r0 / dr + 0.5);
	if (N < 1) {
		N = 1;
	}
	dr = r0 / 32.0;
	r = 0.0;
	do {
		if (r + dr > r0) {
			dr = r0 - r;
			done = 1;
		}
		y0 = y;
		z0 = z;
		dy1 = fy(y, z, r) * dr;
		dz1 = fz(y, z, r) * dr;
		y += 0.5 * dy1;
		z += 0.5 * dz1;
		if (y <= 0.0) {
			y = 0.0;
			break;
		}
		dy2 = fy(y, z, r + 0.5 * dr) * dr;
		dz2 = fz(y, z, r + 0.5 * dr) * dr;
		y = y0 + 0.5 * dy2;
		z = z0 + 0.5 * dz2;
		if (y <= 0.0) {
			y = 0.0;
			break;
		}
		dy3 = fy(y, z, r + 0.5 * dr) * dr;
		dz3 = fz(y, z, r + 0.5 * dr) * dr;
		y = y0 + dy3;
		z = z0 + dz3;
		if (y <= 0.0) {
			y = 0.0;
			break;
		}
		dy4 = fy(y, z, r + dr) * dr;
		dz4 = fz(y, z, r + dr) * dr;
		y = y0 + (dy1 + dy4 + 2.0 * (dy3 + dy2)) / 6.0;
		z = z0 + (dz1 + dz4 + 2.0 * (dz3 + dz2)) / 6.0;
		if (y <= 0.0) {
			y = 0.0;
			break;
		}
		r += dr;
	} while (done == 0);
	if (y < 0.0) {
		return 0.0;
	} else {
		return y;
	}
}
