/*
 * hydro.hpp
 *
 *  Created on: Feb 3, 2015
 *      Author: dmarce1
 */

#ifndef HYDRO_HPP_
#define HYDRO_HPP_

#include "defs.hpp"

class hydro_vars {
public:
	static constexpr bool positive_field[5] = { true, true, false, false, false };
	static constexpr real cfl_factor = (1.0 / 2.0);
	static constexpr integer N3F = (NX + 1) * (NX + 1) * (NX + 1);
	static constexpr integer nf_hydro = 5;
	static constexpr integer bw = 2;
	static constexpr integer d0_i = 0;
	static constexpr integer et_i = 1;
	static constexpr integer s0_i = 2;
	static constexpr integer sx_i = s0_i + 0;
	static constexpr integer sy_i = s0_i + 1;
	static constexpr integer sz_i = s0_i + 2;
	static constexpr integer ro_i = d0_i;
	static constexpr integer pr_i = et_i;
	static constexpr integer v0_i = s0_i;
	static constexpr integer vx_i = sx_i;
	static constexpr integer vy_i = sy_i;
	static constexpr integer vz_i = sz_i;
	static constexpr real ro_floor = 1.0e-5;
	static constexpr real fgamma = 7.0 / 4.0;
private:
	std::vector<std::vector<real>> dU;
	std::vector<std::vector<real>> U0;
	std::vector<std::vector<real>> U;
	std::array<std::vector<std::vector<real>>, NDIM> dUx;
	std::array<std::vector<std::vector<real>>, NDIM> Ux0;
	std::array<std::vector<std::vector<real>>, NDIM> Ux;
	std::array<std::vector<std::vector<real>>, NDIM> FaceFlux;
	std::array<std::vector<std::vector<real>>, NDIM> CenterFlux;
	std::array<std::vector<real>, NDIM> Xc;
	std::array<std::vector<real>, NDIM> Xf;
	std::vector<real> r;
	real dx, amax;

	static integer ind3d(integer j, integer k, integer l, integer stride = NX);
public:
	real cell_mass(integer) const;
	hydro_vars();
	void dump_data(std::vector<double>::iterator& l, integer index) const;
	void initialize(real xcorner, real ycorner, real zcorner, real dx);
	void get_boundary(std::vector<real>::iterator i, integer d) const;
	void set_boundary(integer d, std::vector<real>::iterator*iptr = nullptr);
	void store();
	real compute_du();
	void update(real, integer, const std::vector<real>& phi, const std::vector<real>& gx, const std::vector<real>& gy,
			const std::vector<real>& gz);
	bool needs_refinement() const;
	void restrict_pack(std::vector<real>&, std::vector<bool> amr_dirs) const;
	void restrict_unpack(std::vector<real>&, integer ci, std::vector<bool> amr_directions);
	void prolong_pack(integer dir, integer child, std::vector<real>::iterator i) const;
	void prolong_unpack(integer dir, std::vector<real>::iterator);

};

#endif /* HYDRO_HPP_ */
