/*
 * node_server.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "node_server.hpp"
#include "silo_output.hpp"
#include <mutex>
#include "key.hpp"
#include "exafmm.hpp"

exafmm_kernel exafmm;

constexpr integer WAITING = 0;
constexpr integer READY = 1;
constexpr integer COMPLETE = 2;
constexpr integer lb0 = 0;
constexpr integer lb1 = 0;
constexpr integer lb2 = FMM_NX - FMM_BW;
constexpr integer ub0 = FMM_BW - 1;
constexpr integer ub1 = FMM_NX - 1;
constexpr integer ub2 = FMM_NX - 1;
constexpr integer center_dir = 13;

hpx::id_type node_server::output;

constexpr std::array<bool, NNEIGHBOR> is_face = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };
constexpr std::array<integer, NNEIGHBOR> which_face = { -1, -1, -1, -1, +4, -1, -1, -1, -1,

-1, +2, -1, +0, -1, +1, -1, +3, -1,

-1, -1, -1, -1, +5, -1, -1, -1, -1 };
constexpr std::array<bool, NNEIGHBOR> is_edge = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0 };
constexpr std::array<bool, NNEIGHBOR> is_vert = { 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1 };

constexpr std::array<integer, NNEIGHBOR> cbnd_size = { 1, 2, 1, 2, 4, 2, 1, 2, 1, 2, 4, 2, 4, 8, 4, 2, 4, 2, 1, 2, 1, 2, 4, 2, 1, 2, 1

};

constexpr std::array<integer, NNEIGHBOR> z_off = { -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		+FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW };

constexpr std::array<integer, NNEIGHBOR> y_off = { -FMM_BW, -FMM_BW, -FMM_BW, 0, 0, 0, +FMM_BW, +FMM_BW, +FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, 0, 0, 0, +FMM_BW,
		+FMM_BW, +FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, 0, 0, 0, +FMM_BW, +FMM_BW, +FMM_BW };

constexpr std::array<integer, NNEIGHBOR> x_off = { -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW,
		-FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW };

constexpr std::array<integer, NNEIGHBOR> zlb = { lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb1, lb1, lb1, lb1, lb1, lb1, lb1, lb1, lb1, lb2, lb2, lb2, lb2,
		lb2, lb2, lb2, lb2, lb2 };

constexpr std::array<integer, NNEIGHBOR> ylb = { lb0, lb0, lb0, lb1, lb1, lb1, lb2, lb2, lb2, lb0, lb0, lb0, lb1, lb1, lb1, lb2, lb2, lb2, lb0, lb0, lb0, lb1,
		lb1, lb1, lb2, lb2, lb2 };

constexpr std::array<integer, NNEIGHBOR> xlb = { lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0,
		lb1, lb2, lb0, lb1, lb2 };

constexpr std::array<integer, NNEIGHBOR> zub = { ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub1, ub1, ub1, ub1, ub1, ub1, ub1, ub1, ub1, ub2, ub2, ub2, ub2,
		ub2, ub2, ub2, ub2, ub2 };

constexpr std::array<integer, NNEIGHBOR> yub = { ub0, ub0, ub0, ub1, ub1, ub1, ub2, ub2, ub2, ub0, ub0, ub0, ub1, ub1, ub1, ub2, ub2, ub2, ub0, ub0, ub0, ub1,
		ub1, ub1, ub2, ub2, ub2 };

constexpr std::array<integer, NNEIGHBOR> xub = { ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0,
		ub1, ub2, ub0, ub1, ub2 };

constexpr std::array<integer, NNEIGHBOR> bnd_size = { FMM_BW * FMM_BW * FMM_BW, FMM_BW * FMM_BW * FMM_NX, FMM_BW * FMM_BW * FMM_BW, FMM_BW * FMM_NX * FMM_BW,
		FMM_BW * FMM_NX * FMM_NX, FMM_BW * FMM_NX * FMM_BW, FMM_BW * FMM_BW * FMM_BW, FMM_BW * FMM_BW * FMM_NX, FMM_BW * FMM_BW * FMM_BW, FMM_NX * FMM_BW
				* FMM_BW, FMM_NX * FMM_BW * FMM_NX, FMM_NX * FMM_BW * FMM_BW, FMM_NX * FMM_NX * FMM_BW, FMM_NX * FMM_NX * FMM_NX, FMM_NX * FMM_NX * FMM_BW,
		FMM_NX * FMM_BW * FMM_BW, FMM_NX * FMM_BW * FMM_NX, FMM_NX * FMM_BW * FMM_BW, FMM_BW * FMM_BW * FMM_BW, FMM_BW * FMM_BW * FMM_NX, FMM_BW * FMM_BW
				* FMM_BW, FMM_BW * FMM_NX * FMM_BW, FMM_BW * FMM_NX * FMM_NX, FMM_BW * FMM_NX * FMM_BW, FMM_BW * FMM_BW * FMM_BW, FMM_BW * FMM_BW * FMM_NX,
		FMM_BW * FMM_BW * FMM_BW };

constexpr integer octlb0 = 0;
constexpr integer octlb1 = FMM_NX / 2;
constexpr integer octub0 = FMM_NX / 2 - 1;
constexpr integer octub1 = FMM_NX - 1;

constexpr std::array<integer, NCHILD> oct_zlb = { octlb0, octlb0, octlb0, octlb0, octlb1, octlb1, octlb1, octlb1 };

constexpr std::array<integer, NCHILD> oct_ylb = { octlb0, octlb0, octlb1, octlb1, octlb0, octlb0, octlb1, octlb1 };

constexpr std::array<integer, NCHILD> oct_xlb = { octlb0, octlb1, octlb0, octlb1, octlb0, octlb1, octlb0, octlb1 };

constexpr std::array<integer, NCHILD> oct_zub = { octub0, octub0, octub0, octub0, octub1, octub1, octub1, octub1 };

constexpr std::array<integer, NCHILD> oct_yub = { octub0, octub0, octub1, octub1, octub0, octub0, octub1, octub1 };

constexpr std::array<integer, NCHILD> oct_xub = { octub0, octub1, octub0, octub1, octub0, octub1, octub0, octub1 };

integer ind4d(integer p, integer j, integer k, integer l, integer stride = FMM_NX) {
	return l + stride * (k + stride * (j + stride * p));
}

bool location_is_phys_bnd(integer level, const std::array<integer, NDIM>& loc, integer d) {
	if (dir_x[d] == +1 && loc[0] == (1 << level) - 1) {
		return true;
	} else if (dir_x[d] == -1 && loc[0] == 0) {
		return true;
	} else if (dir_y[d] == +1 && loc[1] == (1 << level) - 1) {
		return true;
	} else if (dir_y[d] == -1 && loc[1] == 0) {
		return true;
	} else if (dir_z[d] == +1 && loc[2] == (1 << level) - 1) {
		return true;
	} else if (dir_z[d] == -1 && loc[2] == 0) {
		return true;
	} else {
		return false;
	}
}

bool node_server::is_phys_bound(integer dir) const {
	return location_is_phys_bnd(level, location, dir);
}

hpx::future<std::vector<real>> node_server::get_expansions(integer c) const {
	return hpx::async([=]() {
		std::vector<real> lc(FMM_PP * FMM_N3 / NCHILD);
		for (integer p = 0; p != FMM_PP; ++p) {
			for (integer j = oct_xlb[c]; j <= oct_xub[c]; ++j) {
				for (integer k = oct_ylb[c]; k <= oct_yub[c]; ++k) {
					for (integer l = oct_zlb[c]; l <= oct_zub[c]; ++l) {
						lc[ind4d(p, j - oct_xlb[c], k - oct_ylb[c], l - oct_zlb[c], FMM_NX / 2)] =
						L[ind4d(p,j,k,l)];
					}
				}
			}
		}
		return std::move(lc);
	});
}

hpx::future<std::vector<real>> node_server::get_multipoles() {
	if (is_leaf && FMM_P != 0) {
		M = hydro_vars->get_gravity_sources(this_rk);
	}
	return hpx::async(hpx::launch::deferred, [=]() {
		std::vector<real> m_out(FMM_PP * FMM_N3 / NCHILD, 0.0);
		std::vector<real> m_in(FMM_PP * FMM_N3 / NCHILD);
		for( integer c = 0; c != NCHILD; ++c) {
			std::array<real,NDIM> dist;

			for (integer p = 0; p != FMM_PP; ++p) {
				for (integer j = 0; j != FMM_NX; j += 2) {
					for (integer k = 0; k != FMM_NX; k += 2) {
						for (integer l = 0; l!= FMM_NX; l += 2) {
							auto tmp = M[ind4d(p,j+(c&1),k+((c>>1)&1),l+((c>>2)&1))];
							auto i = ind4d(p, j>>1, k>>1, l>>1, FMM_NX>>1);
							m_in[i] = tmp;
						}
					}
				}
			}
			for( integer d = 0; d != NDIM; ++d ) {
				if( (c >> d) & 1) {
					dist[d] = -0.5*dx;
				} else {
					dist[d] = +0.5*dx;
				}
			}
			constexpr auto sz = FMM_N3 / NCHILD;
			exafmm.M2M(m_out, m_in, dist, sz);
		}
		return std::move(m_out);
	});
}

bool node_server::refine_me() const {
	if (level < MAXLEVEL) {
		return hydro_vars->refinement_needed();
	} else {
		return false;
	}
}

void node_server::derefine(bool self) {
	hpx::future<void> fut;
	std::vector<hpx::future<void>> cfut(NCHILD);

	if (!is_leaf) {
		std::vector<hpx::future<hpx::id_type>> id_test(NCHILD);
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cfut[ci] = child_id[ci].destroy();
		}
		auto fut = hpx::when_all(cfut).then([](hpx::future<std::vector<hpx::future<void>>>) {});
		fut.get();
		for (integer ci = 0; ci != NCHILD; ++ci) {
			std::array<integer, NDIM> cloc(location);
			for (integer d = 0; d != NDIM; ++d) {
				cloc[d] <<= 1;
				cloc[d] += (ci >> d) & integer(1);
			}
			id_test[ci] = hpx::unregister_id_with_basename("fmmx_node", location_to_key(level + 1, cloc));
		}
		for (integer ci = 0; ci != NCHILD; ++ci) {
			if (id_test[ci].get() != child_id[ci]) {
				printf("Failed to unregister id_with_basename\n");
				abort();
			}
			child_id[ci] = hpx::invalid_id;
		}
		is_leaf = true;
	}
	if (self) {
		my_id = hpx::invalid_id;
		for (int i = 0; i != NNEIGHBOR; ++i) {
			neighbor_id[i] = hpx::invalid_id;
		}
	}
}

void node_server::set_me(hpx::id_type me) {
	my_id = std::move(me);
}

void node_server::refine_proper() {
	std::vector<hpx::future<hpx::id_type>> cids(NCHILD);
	std::vector<hpx::future<void>> cfut(NCHILD);
	{
		boost::lock_guard<hpx::lcos::local::spinlock> lock(R_lock);
		if (!is_leaf) {
			return;
		}
		printf("Refining tree at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);
		for (integer ci = 0; ci != NCHILD; ++ci) {
			std::array<integer, NDIM> cloc(location);
			for (integer d = 0; d != NDIM; ++d) {
				cloc[d] <<= 1;
				cloc[d] += (ci >> d) & integer(1);
			}
			auto locality = key_to_locality(location_to_key(level + 1, cloc));
			cids[ci] = hpx::new_<node_server>(locality, my_id, level + 1, cloc);
		}
		is_leaf = false;
	}
	hpx::wait_all(cids);
	std::vector<hpx::future<bool>> id_test(NCHILD);
	for (integer ci = 0; ci != NCHILD; ++ci) {
		std::array<integer, NDIM> cloc(location);
		for (integer d = 0; d != NDIM; ++d) {
			cloc[d] <<= 1;
			cloc[d] += (ci >> d) & integer(1);
		}
		child_id[ci] = cids[ci].get();
		id_test[ci] = hpx::register_id_with_basename("fmmx_node", child_id[ci], location_to_key(level + 1, cloc));
	}
	for (integer ci = 0; ci != NCHILD; ++ci) {
		if (!id_test[ci].get()) {
			printf("Failed to register id_with_basename\n");
			abort();
		}
	}
	for (integer ci = 0; ci != NCHILD; ++ci) {
		cfut[ci] = child_id[ci].set_me(child_id[ci]);
	}
	for (integer ci = 0; ci != NCHILD; ++ci) {
		cfut[ci].get();
	}

}

bool node_server::refine() {
	const auto refine_this = refine_me();
	if (is_leaf) {
		if (refine_this) {
			refine_proper();
		}
	} else {
		std::vector<hpx::future<bool>> cfut(NCHILD);
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cfut[ci] = child_id[ci].refine();
		}
		std::vector<bool> refine_neighbor(NNEIGHBOR, false);
		for (integer ci = 0; ci != NCHILD; ++ci) {
			if (cfut[ci].get()) {
				for (integer ni = 0; ni != NNEIGHBOR; ++ni) {
					if (ni == center_dir || neighbor_id[ni] == hpx::invalid_id) {
						continue;
					}
					if ((2 * ((ci >> 0) & 1) - 1 == dir_x[ni]) || (2 * ((ci >> 1) & 1) - 1 == dir_y[ni]) || (2 * ((ci >> 2) & 1) - 1 == dir_z[ni])) {
						refine_neighbor[ni] = true;
					}
				}
			}
		}
		std::vector<hpx::future<void>> nfuts(NNEIGHBOR);
		for (integer ni = 0; ni != NNEIGHBOR; ++ni) {
			if (refine_neighbor[ni]) {
				nfuts[ni] = neighbor_id[ni].refine_proper();
			} else {
				nfuts[ni] = hpx::make_ready_future();
			}
		}
		for (integer ni = 0; ni != NNEIGHBOR; ++ni) {
			nfuts[ni].get();
		}
	}
	return !is_leaf;
}

std::vector<real> node_server::hydro_restrict(integer rk) {
	std::vector<hpx::future<std::vector<real>>>cfuts(NCHILD);
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cfuts[ci] = child_id[ci].hydro_restrict(rk);
		}
		for (integer ci = 0; ci != NCHILD; ++ci) {
			auto data = cfuts[ci].get();
			hydro_vars->restrict_unpack(rk, data, ci);
		}
	}
	return hydro_vars->restrict_pack(rk);
}

void node_server::hydro_amr_prolong(integer rk) {

	for (integer f = 0; f != 2 * NDIM; ++f) {
		hydro_vars->set_amr(f, false);
		hydro_vars->set_child_amr(f, false);
	}
	for (integer dir = 0; dir != NNEIGHBOR; ++dir) {
		if (is_face[dir] && !is_phys_bound(dir)) {
			if (neighbor_id[dir] == hpx::invalid_id) {
				hydro_vars->set_amr(which_face[dir], true);
			}
			if (!is_leaf && neighbor_is_leaf[dir]) {
				hydro_vars->set_child_amr(which_face[dir], true);
			}
		}
	}
	hydro_vars->set_compute(is_leaf);
	std::vector<hpx::future<void>> child_exe(NCHILD);
	if (!is_leaf) {
		for (integer i = 0; i != NCHILD; ++i) {
			child_exe[i] = child_id[i].hydro_amr_prolong(rk);
		}
	}
	std::vector<hpx::future<std::vector<real>>>bnd_futs(NNEIGHBOR);
	for (integer dir = 0; dir != NNEIGHBOR; ++dir) {
		bnd_futs[dir] = hpx::make_ready_future(std::vector<real>(0));
		if (dir != center_dir) {
			if (neighbor_id[dir] == hpx::invalid_id && !is_phys_bound(dir) && is_face[dir]) {
				integer ci = 0;
				ci |= ((location[0] & 1) << 0);
				ci |= ((location[1] & 1) << 1);
				ci |= ((location[2] & 1) << 2);
				bnd_futs[dir] = parent_id.hydro_get_amr_bnd(rk, dir, ci);
			}
		}
	}
	for (integer dir = 0; dir != NNEIGHBOR; ++dir) {
		auto data = bnd_futs[dir].get();
		if (data.size()) {
			hydro_vars->unpack_amr_boundary(rk, data, dir);
		}
	}
	if (!is_leaf) {
		for (integer i = 0; i != NCHILD; ++i) {
			child_exe[i].get();
		}
	}
}

void node_server::hydro_exchange(integer rk, exchange_type type) {
	std::vector<hpx::future<void>> child_exe(NCHILD);
	if (!is_leaf) {
		for (integer i = 0; i != NCHILD; ++i) {
			child_exe[i] = child_id[i].hydro_exchange(rk,type);
		}
	}
	std::vector<hpx::future<std::vector<real>>>bnd_futs(NNEIGHBOR);
	for (integer dir = 0; dir != NNEIGHBOR; ++dir) {
		bnd_futs[dir] = hpx::make_ready_future(std::vector<real>(0));
		if (dir != center_dir) {
			if (neighbor_id[dir] != hpx::invalid_id) {
				bnd_futs[dir] = neighbor_id[dir].hydro_get_bnd(rk, NNEIGHBOR - dir - 1, type);
			}
		}
	}
	for (integer dir = 0; dir != NNEIGHBOR; ++dir) {
		auto data = bnd_futs[dir].get();
		if (data.size()) {
			if (neighbor_id[dir] != hpx::invalid_id) {
				switch (type) {
				case HYDRO_BND:
					hydro_vars->unpack_boundary(rk, data, dir);
					break;
				case GRAV_BND:
					hydro_vars->unpack_phi_boundary(rk, data, dir);
					break;
				}
			}
		}
	}
	for (integer dir = 0; dir != NNEIGHBOR; ++dir) {
		if (is_face[dir] && neighbor_id[dir] == hpx::invalid_id && is_phys_bound(dir)) {
			switch (type) {
			case HYDRO_BND:
				hydro_vars->enforce_physical_boundaries(rk, which_face[dir]);
				break;
			case GRAV_BND:
				hydro_vars->enforce_physical_phi_boundaries(rk, which_face[dir]);
				break;
			}
		}
	}
	if (!is_leaf) {
		for (integer i = 0; i != NCHILD; ++i) {
			child_exe[i].get();
		}
	}
}

std::vector<real> node_server::hydro_get_bnd(integer rk, integer dir, exchange_type type) {
	std::vector<real> v;
	switch (type) {
	case HYDRO_BND:
		v= hydro_vars->pack_boundary(rk, dir);
		break;
	case GRAV_BND:
		v= hydro_vars->pack_phi_boundary(rk, dir);
		break;
	}
	return v;
}

std::vector<real> node_server::hydro_get_amr_bnd(integer rk, integer dir, integer ci) {
	return hydro_vars->pack_amr_boundary(rk, dir, ci);
}

hpx::future<std::vector<real>> node_server::get_fmm_boundary(integer d) const {
	return hpx::async(hpx::launch::deferred, [=]() {
		std::vector<real> bnd;
		bnd.resize(FMM_PP * bnd_size[d]);
		auto i = bnd.begin();
		for (integer p = 0; p != FMM_PP; ++p) {
			for (integer j = xlb[d]; j <= xub[d]; ++j) {
				for (integer k = ylb[d]; k <= yub[d]; ++k) {
					for (integer l = zlb[d]; l <= zub[d]; ++l) {
						*i = M[ind4d(p, j, k, l)];
						++i;
					}
				}
			}
		}
		return bnd;
	});
}

void node_server::set_fmm_boundary(hpx::future<std::vector<real>> f, integer d) {
	auto data = f.get();
	M2L(data, d, bnd_size[d]);
	if (++fmm_neighbor_done_cnt == NNEIGHBOR + 1) {
		send_expansions();
	}
}

void node_server::set_expansions(hpx::future<std::vector<real>> f) {
	L2L(f.get());
	if (++fmm_neighbor_done_cnt == NNEIGHBOR + 1) {
		send_expansions();
	}
}

void node_server::set_multipoles(hpx::future<std::vector<real>> f, integer ci) {
	M2M(f.get(), ci);
	if (++fmm_child_done_cnt == NCHILD) {
		send_multipoles();
	}
//	integer cnt = neighbor_done_cnt;
//	if (neighbor_done_cnt == NNEIGHBOR + 1) {
//		send_expansions();
//	}
}

void node_server::M2M(const std::vector<real>& mptr, integer c) {
	auto i = mptr.begin();
	for (integer p = 0; p != FMM_PP; ++p) {
		for (integer j = oct_xlb[c]; j <= oct_xub[c]; ++j) {
			for (integer k = oct_ylb[c]; k <= oct_yub[c]; ++k) {
				for (integer l = oct_zlb[c]; l <= oct_zub[c]; ++l) {
					M[ind4d(p, j, k, l)] = *i;
					++i;
				}
			}
		}
	}
}

void node_server::L2L(const std::vector<real>& l_in) {
	std::vector<real> l_out(FMM_PP * FMM_N3 / NCHILD);
	for (integer c = 0; c != NCHILD; ++c) {
		std::array<real, NDIM> dist;
		for (integer d = 0; d != NDIM; ++d) {
			if ((c >> d) & 1) {
				dist[d] = +0.5 * dx;
			} else {
				dist[d] = -0.5 * dx;
			}
		}
		constexpr auto sz = FMM_N3 / NCHILD;
		exafmm.L2L(l_out, l_in, dist, sz);
		for (integer p = 0; p != FMM_PP; ++p) {
			for (integer j = 0; j != FMM_NX; j += 2) {
				for (integer k = 0; k != FMM_NX; k += 2) {
					for (integer l = 0; l != FMM_NX; l += 2) {
						L[ind4d(p, j + (c & 1), k + ((c >> 1) & 1), l + ((c >> 2) & 1))] += l_out[ind4d(p, j >> 1, k >> 1, l >> 1, FMM_NX / 2)];
					}
				}
			}
		}

	}
}

void node_server::send_multipoles() {
	auto mul_fut = get_multipoles();
	integer j = 0;
	for (integer i = 0; i != NDIM; ++i) {
		j |= ((location[i] & 1) << i);
	}
	parent_id.set_multipoles(mul_fut, j);
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (d == center_dir) {
			continue;
		}
		auto bnd_fut = get_fmm_boundary(d);
		neighbor_id[d].set_fmm_boundary(bnd_fut, NNEIGHBOR - 1 - d);
	}
	M2L(M, center_dir, FMM_N3);
	if (++fmm_neighbor_done_cnt == NNEIGHBOR + 1) {
		send_expansions();
	}
	fmm_child_done_cnt = 0;

}

void node_server::hydro_next_u(integer rk, real dt) {
	std::vector<hpx::future<void>> cfuts(NCHILD);
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cfuts[ci] = child_id[ci].hydro_next_u(rk, dt);
		}
	}
	hydro_vars->next_u(rk, dt);
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cfuts[ci].get();
		}
	}

}

std::pair<real, std::vector<real>> node_server::hydro_next_du(integer rk) {
	std::pair<real, std::vector<real>> rc;
	std::vector<hpx::future<std::pair<real, std::vector<real>>> >cfuts(NCHILD);
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cfuts[ci] = child_id[ci].hydro_next_du(rk);
		}
	}
	rc.first = 1.0e+99;
	std::vector<std::vector<real>> cflux(NCHILD);
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			auto tmp = cfuts[ci].get();
			cflux[ci] = std::move(tmp.second);
			rc.first = std::min(rc.first, tmp.first);
		}
	}
	rc.second = hydro_vars->flux_correct_pack(rk);
	auto this_dt = hydro_vars->next_du(rk, cflux);
	rc.first = std::min(rc.first, this_dt);

	return rc;

}

void node_server::hydro_project(integer rk) {
	std::vector<hpx::future<void>> cfuts(NCHILD);
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cfuts[ci] = child_id[ci].hydro_project(rk);
		}
	}
	hydro_vars->apply_limiter(rk);
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cfuts[ci].get();
		}
	}

}

void node_server::send_expansions() {
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			auto f = get_expansions(ci);
			child_id[ci].set_expansions(f);
		}
	}
	exe_promise->set_value();
}

void node_server::M2L(const std::vector<real>& m, integer d, integer N) {
//printf( "------------%li %li %li %li - %li %li %li\n", level,location[0], location[1], location[2], d % 3, (d/3)%3, d/9);
	if (FMM_P == 0) {
		return;
	}
	auto ub = [=](integer ia) {
		if( level != 0) {
			ia = ((ia+FMM_NX)/2+1)*2 - FMM_NX + 1;
			return std::min(ia,FMM_NX-integer(1));
		} else {
			return FMM_NX - integer(1);
		}
	};

	auto lb = [=](integer ia) {
		if( level != 0 ) {
			ia = ((ia+FMM_NX)/2-1)*2 - FMM_NX;
			return std::max(ia,integer(0));
		} else {
			return integer(0);
		}
	};
#ifdef DIRECT_GRAVITY
	const integer max_d0 =1;
#else
	const integer max_d0 = (is_leaf || neighbor_is_leaf[d]) ? 0 : 1;
#endif
	integer is = 0;
	std::vector<integer> list(level == 0 ? FMM_N3 : 216);
	std::vector<real> l_out;
	std::vector<real> m_in(FMM_PP);
	std::array<std::vector<real>, NDIM> dist;
	for (integer i = 0; i != NDIM; ++i) {
		dist[i].resize(level == 0 ? FMM_N3 : 216);
	}
	integer max_list_size = level == 0 ? FMM_N3 : 216;

	max_list_size = ((max_list_size - 1) / 64 + 1) * 64;
	std::vector<real> L_r(max_list_size), L_i(max_list_size);
	std::vector<real> Ynm(FMM_PP * max_list_size);

	for (integer j1 = xlb[d] + x_off[d]; j1 <= xub[d] + x_off[d]; ++j1) {
		const integer jlb = lb(j1);
		const integer jub = ub(j1);
		for (integer k1 = ylb[d] + y_off[d]; k1 <= yub[d] + y_off[d]; ++k1) {
			const integer klb = lb(k1);
			const integer kub = ub(k1);
			for (integer l1 = zlb[d] + z_off[d]; l1 <= zub[d] + z_off[d]; ++l1) {
				const integer llb = lb(l1);
				const integer lub = ub(l1);
				integer cnt = 0;
				auto i_list = list.begin();
				auto i_xdist = dist[0].begin();
				auto i_ydist = dist[1].begin();
				auto i_zdist = dist[2].begin();

				for (integer j0 = jlb; j0 <= jub; ++j0) {
					for (integer k0 = klb; k0 <= kub; ++k0) {
						const integer tmp = std::max(std::abs(j1 - j0), std::abs(k1 - k0));
						for (integer l0 = llb; l0 <= lub; ++l0) {
							const integer max_d = std::max(tmp, std::abs(l1 - l0));
							if (max_d > max_d0) {
								*i_list++ = ind4d(0, j0, k0, l0);
								*i_xdist++ = real(j0 - j1) * dx;
								*i_ydist++ = real(k0 - k1) * dx;
								*i_zdist++ = real(l0 - l1) * dx;
								++cnt;
							}
						}
					}
				}
				auto i_m = m.begin() + is;
				for (integer p = 0; p < FMM_PP - 1; ++p) {
					m_in[p] = *i_m;
					i_m += N;
				}
				m_in[FMM_PP - 1] = *i_m;
				l_out.resize(FMM_PP * cnt);
				exafmm.M2L(l_out, m_in, dist, cnt, L_r, L_i, Ynm);
				{
					boost::unique_lock<hpx::lcos::local::spinlock> tmp(L_lock);
					for (integer p = 0; p != FMM_PP; ++p) {
						for (integer i = 0; i != cnt; ++i) {
							L[p * FMM_N3 + list[i]] += l_out[p * cnt + i];
						}
					}
				}
				++is;

			}
		}
	}
//	assert(is == N);
}

std::vector<node_client> node_server::get_children() const {
	std::vector<node_client> rc(NCHILD);
	std::copy(child_id.begin(), child_id.end(), rc.begin());
	return rc;
}

void node_server::get_tree(std::vector<node_client> my_neighbors) {
	for (integer i = 0; i != NNEIGHBOR; ++i) {
		neighbor_id[i] = my_neighbors[i];
	}
	neighbor_id[center_dir] = my_id;
	integer cnt = 0;
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (neighbor_id[d] != hpx::invalid_id) {
			++cnt;
		}
	}
	printf("Getting tree at grid %li - %li %li %li - %li\n", level, location[2], location[1], location[0], cnt);
	std::vector<hpx::future<std::vector<node_client>>>futs(NNEIGHBOR);
	std::vector<node_client> neighborhood(NNEIGHBOR * NCHILD);
	for (integer dir = 0; dir != NNEIGHBOR; dir++) {
		futs[dir] = neighbor_id[dir].get_children();
	}
	for (integer dir = 0; dir != NNEIGHBOR; dir++) {
		auto tmp = futs[dir].get();
		if (tmp[0] != hpx::invalid_id) {
			neighbor_is_leaf[dir] = false;
			for (integer ci = 0; ci != NCHILD; ++ci) {
				const auto xi = 2 * (dir_x[dir] + 1) + ((ci >> 0) & 1);
				const auto yi = 2 * (dir_y[dir] + 1) + ((ci >> 1) & 1);
				const auto zi = 2 * (dir_z[dir] + 1) + ((ci >> 2) & 1);
				neighborhood[xi * 36 + yi * 6 + zi] = tmp[ci];
			}
		} else {
			neighbor_is_leaf[dir] = true;
		}
	}
	if (!is_leaf) {
		std::vector<hpx::future<void>> child_futs(NCHILD);
		for (integer ci = 0; ci != NCHILD; ci++) {
			std::vector<node_client> child_neighbors(NNEIGHBOR);
			for (integer dir = 0; dir != NNEIGHBOR; ++dir) {
				const auto xi = (dir_x[dir] + 2) + ((ci >> 0) & 1);
				const auto yi = (dir_y[dir] + 2) + ((ci >> 1) & 1);
				const auto zi = (dir_z[dir] + 2) + ((ci >> 2) & 1);
				child_neighbors[dir] = neighborhood[xi * 36 + yi * 6 + zi];
			}
			child_futs[ci] = child_id[ci].get_tree(child_neighbors);
		}
		hpx::wait_all(child_futs);
	}
	reset();
}

void node_server::init_t0() {
	real cx, cy, cz, dx;
	dx = 1.0 / real(std::pow(2, level) * FMM_NX);
	cx = real(location[0]) / real(std::pow(2, level));
	cy = real(location[1]) / real(std::pow(2, level));
	cz = real(location[2]) / real(std::pow(2, level));
	for (integer p = 1; p < FMM_PP; ++p) {
		for (integer i = 0; i != FMM_N3; ++i) {
			M[FMM_N3 * p + i] = 0.0;
		}
	}
	if (FMM_P > 0) {
		for (integer j = 0; j != FMM_NX; ++j) {
			real x = cx + (real(j) + .5) * dx - .5;
			for (integer k = 0; k != FMM_NX; ++k) {
				real y = cy + (real(k) + .5) * dx - .5;
				for (integer l = 0; l != FMM_NX; ++l) {
					real z = cz + (real(l) + .5) * dx - .5;
					if (sqrt(x * x + y * y + z * z) < .5) {
						M[ind4d(0, j, k, l)] = 1.0;
					} else {
						M[ind4d(0, j, k, l)] = 0.0;
					}
				}
			}
		}
	}
}

bool node_server::is_amr(integer dir) const {
	if (!is_face[dir]) {
		return false;
	} else if (neighbor_id[dir] == hpx::invalid_id) {
		return true;
	} else {
		return false;
	}
}

bool node_server::child_is_amr(integer ci, integer dir) const {
	if (!is_face[dir]) {
		return false;
	} else if (2 * ((ci >> 0) & 1) - 1 == -dir_x[dir]) {
		return false;
	} else if (2 * ((ci >> 1) & 1) - 1 == -dir_y[dir]) {
		return false;
	} else if (2 * ((ci >> 2) & 1) - 1 == -dir_z[dir]) {
		return false;
	} else if (neighbor_id[dir] == hpx::invalid_id) {
		return true;
	} else if (neighbor_is_leaf[dir]) {
		return true;
	} else {
		return false;
	}
}

real node_server::execute(integer rk) {
	this_rk = rk;
	std::vector<hpx::future<real>> child_exe(NCHILD);
//	printf("Begin at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);
	if (!is_leaf) {
		for (integer i = 0; i != NCHILD; ++i) {
			child_exe[i] = child_id[i].execute(rk);
		}
	} else {
		for (integer i = 0; i != NCHILD; ++i) {
			child_exe[i] = hpx::make_ready_future(real(1.0e+99));
		}
		send_multipoles();
	}
	wait_all(child_exe);
	exe_promise->get_future().wait();
	assert(fmm_neighbor_done_cnt == NNEIGHBOR + 1);
	hydro_vars->set_gravity_sources(L, this_rk);
	reset();
	return 0.0;
}

node_server::node_server() {
	assert(false);
}

std::vector<double> node_server::get_data() const {
	auto a = hydro_vars->output_data();
	std::vector<double> b(a.size());
	for (std::size_t i = 0; i != a.size(); ++i) {
		b[i] = double(a[i]);
	}
	return b;
}

node_server::node_server(hpx::id_type sid) {

	output = std::move(sid);
	node_client pid = hpx::naming::invalid_id;
	integer lev = 0;
	std::array<integer, NDIM> loc { { 0, 0, 0 } };
	initialize(std::move(pid), std::move(lev), std::move(loc));
}

node_server::node_server(node_client pid, integer lev, std::array<integer, NDIM> loc) {
	initialize(std::move(pid), std::move(lev), std::move(loc));
}

integer node_server::get_node_count() const {
	integer cnt = 1;
	if (!is_leaf) {
		std::vector<hpx::future<integer>> futs(NCHILD);
		for (integer c = 0; c != NCHILD; ++c) {
			futs[c] = child_id[c].get_node_count();
		}
		hpx::wait_all(futs);
		for (integer c = 0; c != NCHILD; ++c) {
			cnt += futs[c].get();
		}
	}
	return cnt;
}

std::list<std::size_t> node_server::get_leaf_list() const {
	std::list<std::size_t> list;
	if (!is_leaf) {
		std::vector<hpx::future<std::list<std::size_t>>>futs(NCHILD);
		for (integer c = 0; c != NCHILD; ++c) {
			futs[c] = child_id[c].get_leaf_list();
		}
		hpx::wait_all (futs);
		for (integer c = 0; c != NCHILD; ++c) {
			list.splice(list.end(), futs[c].get());
		}
	} else {
		list.push_back(key);
	}
	return list;
}

void node_server::reset() {
	if (FMM_P != 0) {
		/*		for (integer i = 0; i != FMM_N3; ++i) {
		 phi[i] = -L[i];
		 gz[i] = -L[i + FMM_N3 * 2];
		 gx[i] = L[i + FMM_N3 * 3] * std::sqrt(real(2));
		 gy[i] = -L[i + FMM_N3 * 1] * std::sqrt(real(2));
		 }*/
	}
	L = std::vector<real>(FMM_PP * FMM_N3, 0.0);
	fmm_neighbor_done_cnt = 0;
	fmm_child_done_cnt = 0;
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (d != center_dir) {
			if (neighbor_id[d] == hpx::invalid_id) {
				++fmm_neighbor_done_cnt;
			}
		}
	}
	if (level == 0) {
		++fmm_neighbor_done_cnt;
	}
	++step_cnt;
	if (step_cnt % 2 == 0) {
		exe_pair.first = hpx::promise<void>();
		exe_promise = &(exe_pair.first);
	} else {
		exe_pair.second = hpx::promise<void>();
		exe_promise = &(exe_pair.second);
	}

}

void node_server::initialize(node_client pid, integer lev, std::array<integer, NDIM> loc) {
	step_cnt = 0;
	//phi = gx = gy = gz = std::vector<real>(FMM_N3);
	M.resize(FMM_PP * FMM_N3, 0.0);
	L.resize(FMM_PP * FMM_N3, 0.0);
	location = loc;
	level = lev;
	is_leaf = true;
	parent_id = pid;
	key = location_to_key(level, std::move(loc));
	dx = 1.0 / real(std::pow(2, level) * FMM_NX);
	real x0 = real(loc[0]) / real(std::pow(2, level));
	real y0 = real(loc[1]) / real(std::pow(2, level));
	real z0 = real(loc[2]) / real(std::pow(2, level));
	hydro_vars = std::make_shared < hydro > (dx, x0, y0, z0);
	hydro_vars->initialize();
	reset();
}

node_server::~node_server() {
	if (level == 0) {
		output = hpx::invalid_id;
	}
}

