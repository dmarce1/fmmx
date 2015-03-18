/*
 * math.cpp
 *
 *  Created on: Mar 7, 2015
 *      Author: dmarce1
 */

#include <cassert>
#include "math.hpp"


real LegendreP(integer n, real x) {
	switch (n) {
	case 0:
		return real(1);
	case 1:
		return x;
	case 2:
		return (real(3) * x * x - real(1)) / real(2);
	case 3:
		return (real(5) * x * x - real(3)) * (x / real(2));
	case 4:
		assert(false);
	}
}

real dLegendreP_dx(integer n, real x) {
	switch (n) {
	case 0:
		return real(0);
	case 1:
		return real(1);
	case 2:
		return real(3) * x;
	case 3:
		return (real(15) * x * x - real(3)) / real(2);
	case 4:
		assert(false);
	}
}
