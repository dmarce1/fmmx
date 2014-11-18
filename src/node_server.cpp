/*
 * node_server.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "node_server.hpp"
#include "silo_output.hpp"

#include <hpx/lcos/wait_all.hpp>
#include <mutex>
#include "key.hpp"
#include "exafmm.hpp"

exafmm_kernel exafmm;

constexpr integer WAITING = 0;
constexpr integer READY = 1;
constexpr integer COMPLETE = 2;
constexpr integer lb0 = 0;
constexpr integer lb1 = 0;
constexpr integer lb2 = NX - BW;
constexpr integer ub0 = BW - 1;
constexpr integer ub1 = NX - 1;
constexpr integer ub2 = NX - 1;
constexpr integer center_dir = 13;

hpx::id_type node_server::output;

constexpr std::array<bool, NNEIGHBOR> is_face = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1,
		0, 0, 0, 0 };
constexpr std::array<bool, NNEIGHBOR> is_edge = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0,
		1, 0, 1, 0 };
constexpr std::array<bool, NNEIGHBOR> is_vert = { 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
		0, 1, 0, 1 };

constexpr std::array<integer, NNEIGHBOR> cbnd_size = { 1, 2, 1, 2, 4, 2, 1, 2, 1, 2, 4, 2, 4, 8, 4, 2, 4, 2, 1, 2, 1, 2,
		4, 2, 1, 2, 1

};

constexpr std::array<integer, NNEIGHBOR> z_off = { -BW, -BW, -BW, -BW, -BW, -BW, -BW, -BW, -BW, 0, 0, 0, 0, 0, 0, 0, 0,
		0, +BW, +BW, +BW, +BW, +BW, +BW, +BW, +BW, +BW };

constexpr std::array<integer, NNEIGHBOR> y_off = { -BW, -BW, -BW, 0, 0, 0, +BW, +BW, +BW, -BW, -BW, -BW, 0, 0, 0, +BW,
		+BW, +BW, -BW, -BW, -BW, 0, 0, 0, +BW, +BW, +BW };

constexpr std::array<integer, NNEIGHBOR> x_off = { -BW, 0, +BW, -BW, 0, +BW, -BW, 0, +BW, -BW, 0, +BW, -BW, 0, +BW, -BW,
		0, +BW, -BW, 0, +BW, -BW, 0, +BW, -BW, 0, +BW };

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

constexpr std::array<integer, NNEIGHBOR> bnd_size = { BW * BW * BW, BW * BW * NX, BW * BW * BW,

BW * NX * BW, BW * NX * NX, BW * NX * BW,

BW * BW * BW, BW * BW * NX, BW * BW * BW,

NX * BW * BW, NX * BW * NX, NX * BW * BW,

NX * NX * BW, NX * NX * NX, NX * NX * BW,

NX * BW * BW, NX * BW * NX, NX * BW * BW,

BW * BW * BW, BW * BW * NX, BW * BW * BW,

BW * NX * BW, BW * NX * NX, BW * NX * BW,

BW * BW * BW, BW * BW * NX, BW * BW * BW };

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

void node_server::reset_children() {
	if (is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_status[ci] = WAITING;
		}
	} else {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_status[ci] = COMPLETE;
		}
	}
}

void node_server::reset_parent() {
	parent_status = level == 0 ? COMPLETE : WAITING;
}

void node_server::reset_neighbors() {
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (location_is_phys_bnd(level, location, d)) {
			neighbor_status[d] = COMPLETE;
		} else {
			neighbor_status[d] = WAITING;
		}
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

hpx::future<std::vector<real>> node_server::get_multipoles() const {
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

void node_server::refine() {
	if (my_id == hpx::invalid_id) {
		my_id = hpx::find_id_from_basename("fmmx_node", location_to_key(level, location)).get();
	}
	std::vector<hpx::future<hpx::id_type>> cids(NCHILD);
	for (integer ci = 0; ci != NCHILD; ++ci) {
		std::array<integer, NDIM> cloc(location);
		for (integer d = 0; d != NDIM; ++d) {
			cloc[d] <<= 1;
			cloc[d] += (ci >> d) & integer(1);
		}
		cids[ci] = hpx::new_ < node_server > (hpx::find_here(), get_gid(), level + 1, cloc);

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
}

hpx::future<std::vector<real>> node_server::get_boundary(integer d) const {
	return hpx::async(hpx::launch::deferred, [=]() {
		std::vector<real> bnd(PP * bnd_size[d]);
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
		assert( i == bnd.end() );
		return bnd;
	});
}

void node_server::set_boundary(hpx::future<std::vector<real>> f, integer d) {
	neighbor_futures[d] = std::move(f);
	neighbor_status[d] = READY;
	input_condition.signal();
}

void node_server::set_multipoles(hpx::future<std::vector<real>> f, integer ci) {
	child_futures[ci] = std::move(f);
	child_status[ci] = READY;
	input_condition.signal();
}

void node_server::set_expansions(hpx::future<std::vector<real>> f) {
	parent_future = std::move(f);
	parent_status = READY;
	input_condition.signal();
}

void node_server::wait_for_signal() const {
	input_condition.wait();
}

void node_server::M2M(const std::vector<real>& mptr, integer c) {
	auto i = mptr.begin();
	for (integer p = 0; p != PP; ++p) {
		for (integer j = oct_xlb[c]; j <= oct_xub[c]; ++j) {
			for (integer k = oct_ylb[c]; k <= oct_yub[c]; ++k) {
				for (integer l = oct_zlb[c]; l <= oct_zub[c]; ++l) {
					M[ind4d(p, j, k, l)] = *i;
					L[ind4d(p, j, k, l)] = 0.0;
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
						L[ind4d(p, j + (c & 1), k + ((c >> 1) & 1), l + ((c >> 2) & 1))] += l_in[ind4d(p, j >> 1,
								k >> 1, l >> 1, NX / 2)];
					}
				}
			}
		}

	}
}

void node_server::M2L(const std::vector<real>& m, integer d) {
	//printf( "------------%li %li %li %li - %li %li %li\n", level,location[0], location[1], location[2], d % 3, (d/3)%3, d/9);
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
	const integer N = m.size() / PP;
	integer is = 0;
	std::vector<integer> list(level == 0 ? N3 : 216);
	std::vector<real> l_out;
	std::vector<real> m_in(PP);
	std::array<std::vector<real>, NDIM> dist;
	for (integer i = 0; i != NDIM; ++i) {
		dist[i].resize(level == 0 ? N3 : 216);
	}
	const real delta = 1.0 / real(std::pow(integer(2), level));
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
				exafmm.M2L(l_out, m_in, dist, cnt);
				for (integer p = 0; p != PP; ++p) {
					for (integer i = 0; i != cnt; ++i) {
						L[p * N3 + list[i]] += l_out[p * cnt + i];
					}
				}
				++is;

			}
		}
	}
	assert(is == N);
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
	hpx::wait_all(std::move(f));
	auto j = f.begin();
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (!location_is_phys_bnd(level, location, d)) {
			neighbor_id[d] = j->get();
			++j;
		}
	}
}

void node_server::init_t0() {
	real cx, cy, cz, dx;
	dx = 1.0 / real(std::pow(2, level) * NX);
	cx = real(location[0]) / real(std::pow(2, level));
	cy = real(location[1]) / real(std::pow(2, level));
	cz = real(location[2]) / real(std::pow(2, level));
	for (integer p = 1; p != PP; ++p) {
		for (integer i = 0; i != N3; ++i) {
			M[N3 * p + i] = 0.0;
		}
	}
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

void node_server::execute() {
	bool done;
	printf("Begin at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);
	if (level < MAXLEVEL) {
		refine();
	}
	init_t0();
	get_tree();

	printf("Ascend at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);
	hpx::this_thread::sleep_for(boost::posix_time::milliseconds(1));
	if (!is_leaf) {
		do {
			wait_for_signal();
			done = true;
			for (integer c = 0; c != NCHILD; ++c) {
				switch (child_status[c]) {
				case WAITING:
					done = false;
					break;
				case READY:
					M2M(child_futures[c].get(), c);
					child_status[c] = COMPLETE;
				default:
					break;
				}
			}
		} while (!done);
	}
	auto mul_fut = get_multipoles();
	integer j = 0;
	for (integer i = 0; i != NDIM; ++i) {
		j |= ((location[i] & 1) << i);
	}
	parent_id.set_multipoles(mul_fut, j);
	printf("Exchange at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (d == center_dir) {
			continue;
		}
		auto bnd_fut = get_boundary(d);
		neighbor_id[d].set_boundary(bnd_fut, NNEIGHBOR - 1 - d);
	}
	M2L(M, center_dir);
	if (level > 0) {
		neighbor_status[center_dir] = COMPLETE;
		do {
			wait_for_signal();
			done = true;
			for (integer d = 0; d != NNEIGHBOR; ++d) {
				switch (neighbor_status[d]) {
				case WAITING:
					done = false;
					break;
				case READY:
					M2L(neighbor_futures[d].get(), d);
					neighbor_status[d] = COMPLETE;
				default:
					break;
				}
			}
		} while (!done);
		do {
			wait_for_signal();
			done = true;
			switch (parent_status) {
			case WAITING:
				done = false;
				break;
			case READY:
				L2L(parent_future.get());
				parent_status = COMPLETE;
			default:
				break;
			}
		} while (!done);
	}
	printf("Descend at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);
	if (!is_leaf) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			auto f = get_expansions(ci);
			child_id[ci].set_expansions(f);
		}
	}
	printf("End at grid %li - %li %li %li\n", level, location[2], location[1], location[0]);
	data_ready.signal();
	if (level == 0) {
		std::list<std::size_t> leaf_list = get_leaf_list();
		printf("%li leaves detected by root\n", leaf_list.size());
		auto f = hpx::async<typename silo_output::do_output_action>(output, std::move(leaf_list));
		f.get();
	}
}

node_server::node_server() {
	assert(false);
}

std::vector<real> node_server::get_data() const {
	data_ready.wait();
	printf("Getting data\n");
	std::vector<real> v((2 * PP) * N3);
	auto l = v.begin();
	for (integer i = 0; i != N3; ++i) {
		for (integer p = 0; p != PP; ++p) {
			*l++ = M[N3 * p + i];
		}
		for (integer p = 0; p != PP; ++p) {
			*l++ = L[N3 * p + i];
		}
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
		hpx::wait_all(futs);
		for (integer c = 0; c != NCHILD; ++c) {
			list.splice(list.end(), futs[c].get());
		}
	} else {
		list.push_back(key);
	}
	return list;
}

void node_server::initialize(node_client pid, integer lev, std::array<integer, NDIM> loc) {
	M.resize(PP * N3);
	L.resize(PP * N3);
	neighbors_set = false;
	location = loc;
	level = lev;
	parent_id = node_client(pid);
	is_leaf = true;
	reset_children();
	reset_parent();
	reset_neighbors();
	key = location_to_key(level, std::move(loc));
	dx = 1.0 / real(std::pow(2, level) * NX);
	my_thread = hpx::thread([=]() {
		execute();
	});
}

node_server::~node_server() {
	my_thread.join();
}

