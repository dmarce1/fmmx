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
constexpr integer lb2 = NX - FMM_BW;
constexpr integer ub0 = FMM_BW - 1;
constexpr integer ub1 = NX - 1;
constexpr integer ub2 = NX - 1;
constexpr integer center_dir = 13;

hpx::id_type node_server::output;

constexpr std::array<bool, NNEIGHBOR> is_face = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1,
		0, 0, 0, 0 };
constexpr std::array<integer, NNEIGHBOR> which_face = { -1, -1, -1, -1, +4, -1, -1, -1, -1,

-1, +2, -1, +0, -1, +1, -1, +3, -1,

-1, -1, -1, -1, +5, -1, -1, -1, -1 };
constexpr std::array<bool, NNEIGHBOR> is_edge = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0,
		1, 0, 1, 0 };
constexpr std::array<bool, NNEIGHBOR> is_vert = { 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
		0, 1, 0, 1 };

constexpr std::array<integer, NNEIGHBOR> cbnd_size = { 1, 2, 1, 2, 4, 2, 1, 2, 1, 2, 4, 2, 4, 8, 4, 2, 4, 2, 1, 2, 1, 2,
		4, 2, 1, 2, 1

};

constexpr std::array<integer, NNEIGHBOR> z_off = { -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW,
		-FMM_BW, -FMM_BW, 0, 0, 0, 0, 0, 0, 0, 0, 0, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW, +FMM_BW,
		+FMM_BW, +FMM_BW };

constexpr std::array<integer, NNEIGHBOR> y_off = { -FMM_BW, -FMM_BW, -FMM_BW, 0, 0, 0, +FMM_BW, +FMM_BW, +FMM_BW,
		-FMM_BW, -FMM_BW, -FMM_BW, 0, 0, 0, +FMM_BW, +FMM_BW, +FMM_BW, -FMM_BW, -FMM_BW, -FMM_BW, 0, 0, 0, +FMM_BW,
		+FMM_BW, +FMM_BW };

constexpr std::array<integer, NNEIGHBOR> x_off = { -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW,
		-FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW, -FMM_BW, 0, +FMM_BW,
		-FMM_BW, 0, +FMM_BW };

constexpr std::array<integer, NNEIGHBOR> zlb = { lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb1, lb1, lb1, lb1, lb1,
		lb1, lb1, lb1, lb1, lb2, lb2, lb2, lb2, lb2, lb2, lb2, lb2, lb2 };

constexpr std::array<integer, NNEIGHBOR> ylb = { lb0, lb0, lb0, lb1, lb1, lb1, lb2, lb2, lb2, lb0, lb0, lb0, lb1, lb1,
		lb1, lb2, lb2, lb2, lb0, lb0, lb0, lb1, lb1, lb1, lb2, lb2, lb2 };

constexpr std::array<integer, NNEIGHBOR> xlb = { lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1,
		lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2 };

constexpr std::array<integer, NNEIGHBOR> zub = { ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub1, ub1, ub1, ub1, ub1,
		ub1, ub1, ub1, ub1, ub2, ub2, ub2, ub2, ub2, ub2, ub2, ub2, ub2 };

constexpr std::array<integer, NNEIGHBOR> yub = { ub0, ub0, ub0, ub1, ub1, ub1, ub2, ub2, ub2, ub0, ub0, ub0, ub1, ub1,
		ub1, ub2, ub2, ub2, ub0, ub0, ub0, ub1, ub1, ub1, ub2, ub2, ub2 };

constexpr std::array<integer, NNEIGHBOR> xub = { ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1,
		ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2 };

constexpr std::array<integer, NNEIGHBOR> bnd_size = { FMM_BW * FMM_BW * FMM_BW, FMM_BW * FMM_BW * NX, FMM_BW * FMM_BW
		* FMM_BW, FMM_BW * NX * FMM_BW, FMM_BW * NX * NX, FMM_BW * NX * FMM_BW, FMM_BW * FMM_BW * FMM_BW, FMM_BW
		* FMM_BW * NX, FMM_BW * FMM_BW * FMM_BW, NX * FMM_BW * FMM_BW, NX * FMM_BW * NX, NX * FMM_BW * FMM_BW, NX * NX
		* FMM_BW, NX * NX * NX, NX * NX * FMM_BW, NX * FMM_BW * FMM_BW, NX * FMM_BW * NX, NX * FMM_BW * FMM_BW, FMM_BW
		* FMM_BW * FMM_BW, FMM_BW * FMM_BW * NX, FMM_BW * FMM_BW * FMM_BW, FMM_BW * NX * FMM_BW, FMM_BW * NX * NX,
		FMM_BW * NX * FMM_BW, FMM_BW * FMM_BW * FMM_BW, FMM_BW * FMM_BW * NX, FMM_BW * FMM_BW * FMM_BW };

constexpr integer octlb0 = 0;
constexpr integer octlb1 = NX / 2;
constexpr integer octub0 = NX / 2 - 1;
constexpr integer octub1 = NX - 1;

constexpr std::array<integer, NCHILD> oct_zlb = { octlb0, octlb0, octlb0, octlb0, octlb1, octlb1, octlb1, octlb1 };

constexpr std::array<integer, NCHILD> oct_ylb = { octlb0, octlb0, octlb1, octlb1, octlb0, octlb0, octlb1, octlb1 };

constexpr std::array<integer, NCHILD> oct_xlb = { octlb0, octlb1, octlb0, octlb1, octlb0, octlb1, octlb0, octlb1 };

constexpr std::array<integer, NCHILD> oct_zub = { octub0, octub0, octub0, octub0, octub1, octub1, octub1, octub1 };

constexpr std::array<integer, NCHILD> oct_yub = { octub0, octub0, octub1, octub1, octub0, octub0, octub1, octub1 };

constexpr std::array<integer, NCHILD> oct_xub = { octub0, octub1, octub0, octub1, octub0, octub1, octub0, octub1 };

constexpr std::array<integer, NNEIGHBOR> dir_z = { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
		1, 1, 1, 1, 1, 1, 1 };

constexpr std::array<integer, NNEIGHBOR> dir_y = { -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1,
		-1, 0, 0, 0, 1, 1, 1 };

constexpr std::array<integer, NNEIGHBOR> dir_x = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1,
		-1, 0, 1, -1, 0, 1 };

integer ind4d(integer p, integer j, integer k, integer l, integer stride = NX) {
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

hpx::future<std::vector<real>> node_server::get_expansions(integer c) const {
	return hpx::async([=]() {
		std::vector<real> lc(PP * N3 / NCHILD);
		for (integer p = 0; p != PP; ++p) {
			for (integer j = oct_xlb[c]; j <= oct_xub[c]; ++j) {
				for (integer k = oct_ylb[c]; k <= oct_yub[c]; ++k) {
					for (integer l = oct_zlb[c]; l <= oct_zub[c]; ++l) {
						lc[ind4d(p, j - oct_xlb[c], k - oct_ylb[c], l - oct_zlb[c], NX / 2)] =
						L[ind4d(p,j,k,l)];
					}
				}
			}
		}
		return std::move(lc);
	});
}

hpx::future<std::vector<real>> node_server::get_multipoles() {
	if (is_leaf) {
		M = std::vector<real>(N3 * PP, real(0));
		for (integer i = 0; i != N3; ++i) {
			M[i] = hydro_state.cell_mass(i);
		}
	}
	return hpx::async(hpx::launch::deferred, [=]() {
		std::vector<real> m_out(PP * N3 / NCHILD, 0.0);
		std::vector<real> m_in(PP * N3 / NCHILD);
		for( integer c = 0; c != NCHILD; ++c) {
			std::array<real,NDIM> dist;

			for (integer p = 0; p != PP; ++p) {
				for (integer j = 0; j != NX; j += 2) {
					for (integer k = 0; k != NX; k += 2) {
						for (integer l = 0; l!= NX; l += 2) {
							auto tmp = M[ind4d(p,j+(c&1),k+((c>>1)&1),l+((c>>2)&1))];
							auto i = ind4d(p, j>>1, k>>1, l>>1, NX>>1);
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
			constexpr auto sz = N3 / NCHILD;
			exafmm.M2M(m_out, m_in, dist, sz);
		}
		return std::move(m_out);
	});
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

void node_server::refine() {
	if (my_id == hpx::invalid_id) {
		my_id = hpx::find_id_from_basename("fmmx_node", location_to_key(level, location)).get();
	}
	hpx::future<void> fut;
	std::vector<hpx::future<void>> cfut(NCHILD);
	if (level < MAXLEVEL) {
		std::vector<hpx::future<hpx::id_type>> cids(NCHILD);
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
			cfut[ci] = child_id[ci].refine();
		}
	}
	init_t0();
	if (level < MAXLEVEL) {
		fut = hpx::when_all(cfut).then([](hpx::future<std::vector<hpx::future<void>>>) {});
	} else {
		fut = hpx::make_ready_future();
	}
	fut.get();
	get_tree();
}

hpx::future<std::vector<real>> node_server::get_boundary(integer d) const {
	return hpx::async(hpx::launch::deferred, [=]() {
		std::vector<real> bnd;
		if( is_face[d]) {
			bnd.resize(PP * bnd_size[d] + hydro_vars::nf_hydro * hydro_vars::bw * NX * NX);
		} else {
			bnd.resize(PP * bnd_size[d]);
		}
		auto i = bnd.begin();
		for (integer p = 0; p != PP; ++p) {
			for (integer j = xlb[d]; j <= xub[d]; ++j) {
				for (integer k = ylb[d]; k <= yub[d]; ++k) {
					for (integer l = zlb[d]; l <= zub[d]; ++l) {
						*i = M[ind4d(p, j, k, l)];
						++i;
					}
				}
			}
		}
		if( is_face[d]) {
			hydro_state.get_boundary(bnd.begin() + PP * bnd_size[d], which_face[d]);
		}
		return bnd;
	});
}

void node_server::set_boundary(hpx::future<std::vector<real>> f, integer d) {
	auto data = f.get();
	M2L(data, d, bnd_size[d]);
	if (is_face[d]) {
		if (is_leaf) {
			while (!hydro_updated) {
				hpx::this_thread::yield();
			}
		}
		auto iter = data.begin() + PP * bnd_size[d];
		hydro_state.set_boundary(which_face[d], &iter);
	}
	integer cnt = neighbor_done_cnt;
//	printf("%li %li\n", level, cnt);
	if (++neighbor_done_cnt == NNEIGHBOR + 1) {
		send_expansions();
	}
}

void node_server::set_expansions(hpx::future<std::vector<real>> f) {
	L2L(f.get());
	integer cnt = neighbor_done_cnt;
//	printf("%li %li\n", level, cnt);

	if (++neighbor_done_cnt == NNEIGHBOR + 1) {
		send_expansions();
	}
}

void node_server::set_multipoles(hpx::future<std::vector<real>> f, integer ci) {
	M2M(f.get(), ci);
	if (++child_done_cnt == NCHILD) {
		send_multipoles();
	}
	integer cnt = neighbor_done_cnt;
	if (neighbor_done_cnt == NNEIGHBOR + 1) {
		send_expansions();
	}
}

void node_server::M2M(const std::vector<real>& mptr, integer c) {
	auto i = mptr.begin();
	for (integer p = 0; p != PP; ++p) {
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
	std::vector<real> l_out(PP * N3 / NCHILD);
	for (integer c = 0; c != NCHILD; ++c) {
		std::array<real, NDIM> dist;
		for (integer d = 0; d != NDIM; ++d) {
			if ((c >> d) & 1) {
				dist[d] = -0.5 * dx;
			} else {
				dist[d] = +0.5 * dx;
			}
		}
		constexpr auto sz = N3 / NCHILD;
		exafmm.L2L(l_out, l_in, dist, sz);
		for (integer p = 0; p != PP; ++p) {
			for (integer j = 0; j != NX; j += 2) {
				for (integer k = 0; k != NX; k += 2) {
					for (integer l = 0; l != NX; l += 2) {
						L[ind4d(p, j + (c & 1), k + ((c >> 1) & 1), l + ((c >> 2) & 1))] += l_out[ind4d(p, j >> 1,
								k >> 1, l >> 1, NX / 2)];
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
		auto bnd_fut = get_boundary(d);
		neighbor_id[d].set_boundary(bnd_fut, NNEIGHBOR - 1 - d);
	}
	M2L(M, center_dir, N3);
	neighbor_done_cnt++;
	child_done_cnt = 0;

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
	if (P == 0) {
		return;
	}
	auto ub = [=](integer ia) {
		if( level != 0) {
			ia = ((ia+NX)/2+1)*2 - NX + 1;
			return std::min(ia,NX-integer(1));
		} else {
			return NX - integer(1);
		}
	};

	auto lb = [=](integer ia) {
		if( level != 0 ) {
			ia = ((ia+NX)/2-1)*2 - NX;
			return std::max(ia,integer(0));
		} else {
			return integer(0);
		}
	};

	const integer max_d0 = is_leaf ? 0 : 1;
	integer is = 0;
	std::vector<integer> list(level == 0 ? N3 : 216);
	std::vector<real> l_out;
	std::vector<real> m_in(PP);
	std::array<std::vector<real>, NDIM> dist;
	for (integer i = 0; i != NDIM; ++i) {
		dist[i].resize(level == 0 ? N3 : 216);
	}
	const real delta = 1.0 / real(std::pow(integer(2), level));

	integer max_list_size = level == 0 ? N3 : 216;

	max_list_size = ((max_list_size - 1) / 64 + 1) * 64;
	std::vector<real> L_r(max_list_size), L_i(max_list_size);
	std::vector<real> Ynm(PP * max_list_size);

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
								*i_xdist++ = real(j1 - j0) * dx;
								*i_ydist++ = real(k1 - k0) * dx;
								*i_zdist++ = real(l1 - l0) * dx;
								++cnt;
							}
						}
					}
				}
				auto i_m = m.begin() + is;
				for (integer p = 0; p < PP - 1; ++p) {
					m_in[p] = *i_m;
					i_m += N;
				}
				m_in[PP - 1] = *i_m;
				l_out.resize(PP * cnt);
				exafmm.M2L(l_out, m_in, dist, cnt, L_r, L_i, Ynm);
				{
					boost::unique_lock<hpx::lcos::local::spinlock> tmp(L_lock);
					for (integer p = 0; p != PP; ++p) {
						for (integer i = 0; i != cnt; ++i) {
							L[p * N3 + list[i]] += l_out[p * cnt + i];
						}
					}
				}
				++is;

			}
		}
	}
//	assert(is == N);
}

void node_server::get_tree() {
	std::vector<std::size_t> ids(NNEIGHBOR);
	auto i = ids.begin();
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (!location_is_phys_bnd(level, location, d)) {
			auto this_loc = location;
			this_loc[0] += dir_x[d];
			this_loc[1] += dir_y[d];
			this_loc[2] += dir_z[d];
			*i = location_to_key(level, std::move(this_loc));
			++i;
		}
	}
	auto f = hpx::find_ids_from_basename("fmmx_node", std::vector<std::size_t>(ids.begin(), i));
//	hpx::wait_all(f);
	auto j = f.begin();
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (!location_is_phys_bnd(level, location, d)) {
			neighbor_id[d] = j->get();
			++j;
		}
	}
	printf("Tree at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);

}

void node_server::init_t0() {
	real cx, cy, cz, dx;
	dx = 1.0 / real(std::pow(2, level) * NX);
	cx = real(location[0]) / real(std::pow(2, level));
	cy = real(location[1]) / real(std::pow(2, level));
	cz = real(location[2]) / real(std::pow(2, level));
	for (integer p = 1; p < PP; ++p) {
		for (integer i = 0; i != N3; ++i) {
			M[N3 * p + i] = 0.0;
		}
	}
	if (P > 0) {
		for (integer j = 0; j != NX; ++j) {
			real x = cx + (real(j) + .5) * dx - .5;
			for (integer k = 0; k != NX; ++k) {
				real y = cy + (real(k) + .5) * dx - .5;
				for (integer l = 0; l != NX; ++l) {
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

real node_server::execute(real global_dt, integer rk) {
	bool done, all_poles, poles_sent;
	std::vector<hpx::future<real>> child_exe(NCHILD);
//	printf("Begin at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);
	if (!is_leaf) {
		for (integer i = 0; i != NCHILD; ++i) {
			child_exe[i] = child_id[i].execute(global_dt, rk);
		}
	} else {
		auto fiter = L.begin();
		fiter += N3;
		hydro_state.update(global_dt, rk, phi, gx, gy, gz);
		hydro_updated = true;
		for (integer i = 0; i != NCHILD; ++i) {
			child_exe[i] = hpx::make_ready_future(real(0));
		}
		send_multipoles();
	}
	wait_all(child_exe);
	exe_promise->get_future().wait();
	if (is_leaf) {
		for (integer d = 0; d != NNEIGHBOR; ++d) {
			if (d != center_dir) {
				if (location_is_phys_bnd(level, location, d) && is_face[d]) {
					hydro_state.set_boundary(which_face[d]);
				}
			}
		}
		assert(neighbor_done_cnt == NNEIGHBOR + 1);
		this_dt = hydro_state.compute_du();
	} else {
		if (!is_leaf) {
			this_dt = std::numeric_limits<real>::max();
			for (integer i = 0; i != NCHILD; ++i) {
				this_dt = std::min(this_dt, child_exe[i].get());
			}
		}
	}
	reset();
//	printf("End at grid %li - %li %li %li - dt = %e\n", level, location[2], location[1], location[0], this_dt);
	return this_dt;
}

node_server::node_server() {
	assert(false);
}

std::vector<double> node_server::get_data() const {
//	printf("Getting data\n");
	std::vector<double> v((4 + hydro_vars::nf_hydro) * N3);
	auto l = v.begin();
	for (integer i = 0; i != N3; ++i) {
		*l++ = double(phi[i]);
		*l++ = double(gx[i]);
		*l++ = double(gy[i]);
		*l++ = double(gz[i]);
		hydro_state.dump_data(l, i);
	}
	return std::move(v);
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
	for (integer i = 0; i != N3; ++i) {
		phi[i] = -L[i];
		gz[i] = -L[i + N3 * 2];
		gx[i] = L[i + N3 * 3] * std::sqrt(real(2));
		gy[i] = -L[i + N3 * 1] * std::sqrt(real(2));
	}
	L = std::vector<real>(PP * N3, 0.0);
	hydro_updated = false;
	neighbor_done_cnt = 0;
	child_done_cnt = 0;
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (d != center_dir) {
			if (location_is_phys_bnd(level, location, d)) {
				++neighbor_done_cnt;
			}
		}
	}
	if (level == 0) {
		++neighbor_done_cnt;
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
	phi = gx = gy = gz = std::vector<real>(N3);
	M.resize(PP * N3, 0.0);
	L.resize(PP * N3, 0.0);
	location = loc;
	level = lev;
	is_leaf = true;
	parent_id = pid;
	key = location_to_key(level, std::move(loc));
	dx = 1.0 / real(std::pow(2, level) * NX);
	reset();
	hydro_state.initialize(dx * real(NX * loc[0]), dx * real(NX * loc[1]), dx * real(NX * loc[2]), dx);
}

node_server::~node_server() {
	if (level == 0) {
		output = hpx::invalid_id;
	}
}

