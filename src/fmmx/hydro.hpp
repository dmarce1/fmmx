/*
 * hydro.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: dmarce1
 */

#ifndef HYDRO_HPP_
#define HYDRO_HPP_

#include <valarray>
#include "defs.hpp"
#include "state.hpp"

class hydro {
private:
	real dx, x0, y0, z0;
	std::array<bool, HYDRO_N3> is_amr;
	std::array<bool, HYDRO_N3> is_child_amr;
	std::valarray<std::valarray<std::valarray<std::valarray<real>>> >U;
	std::valarray<std::valarray<std::valarray<std::valarray<real>>>> dU;
	static integer pindex(integer, integer, integer);
	static integer gindex(integer, integer, integer, integer=HYDRO_NX);
	void transform(const std::function<state(const state&)>&, integer rk);
	static state value_at(const std::valarray<state>&, real, real, real);
public:
	std::vector<real> restrict_pack( integer rk ) const;
	void restrict_unpack( integer rk, const std::vector<real>& data, integer ci );
	void set_amr( integer face , bool val=true);
	void set_child_amr( integer face , bool val=true);
 	state get_U_at(real, real, real) const;
	void apply_limiter(integer rk);
	real next_du(integer rk);
	void next_u(integer rk, real dt);
	void initialize();
	void enforce_physical_boundaries(integer rk, integer dir);
	std::vector<real> output_data() const;
	hydro(real dx, real x0, real y0, real z0);
	std::vector<real> pack_boundary(integer rk, integer dir) const;
	void unpack_boundary(integer rk, const std::vector<real>&, integer dir);
	std::vector<real> pack_amr_boundary(integer rk, integer dir, integer ci) const;
	void unpack_amr_boundary(integer rk, const std::vector<real>&, integer dir);
	bool refinement_needed() const;
};

#endif /* HYDRO_HPP_ */
