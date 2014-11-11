/*
 * node_server.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "node_server.hpp"

#include <hpx/lcos/wait_all.hpp>
#include <mutex>

constexpr integer WAITING = 0;
constexpr integer READY = 1;
constexpr integer COMPLETE = 2;
constexpr integer lb0 = 0;
constexpr integer lb1 = 0;
constexpr integer lb2 = NX - BW;
constexpr integer ub0 = BW;
constexpr integer ub1 = NX - 1;
constexpr integer ub2 = NX - 1;
constexpr integer center_dir = 13;

constexpr std::array<integer, NNEIGHBOR> cbnd_size = { 1, 2, 1, 2, 4, 2, 1, 2, 1, 2, 4, 2, 4, 8, 4, 2, 4, 2, 1, 2, 1, 2,
		4, 2, 1, 2, 1

};

constexpr std::array<integer, NNEIGHBOR> x_off = { -BW, -BW, -BW, BW, -BW, -BW, BW, -BW, -BW, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		+BW, +BW, +BW, +BW, +BW, +BW, +BW, +BW, +BW };

constexpr std::array<integer, NNEIGHBOR> y_off = { -BW, -BW, -BW, 0, 0, 0, +BW, +BW, +BW, -BW, -BW, -BW, 0, 0, 0, +BW,
		+BW, +BW, -BW, -BW, -BW, 0, 0, 0, +BW, +BW, +BW };

constexpr std::array<integer, NNEIGHBOR> z_off = { -BW, 0, +BW, -BW, 0, +BW, -BW, 0, +BW, -BW, 0, +BW, -BW, 0, +BW, -BW,
		0, +BW, -BW, 0, +BW, -BW, 0, +BW, -BW, 0, +BW };

constexpr std::array<integer, NNEIGHBOR> xlb = { lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb0, lb1, lb1, lb1, lb1, lb1,
		lb1, lb1, lb1, lb1, lb2, lb2, lb2, lb2, lb2, lb2, lb2, lb2, lb2 };

constexpr std::array<integer, NNEIGHBOR> ylb = { lb0, lb0, lb0, lb1, lb1, lb1, lb2, lb2, lb2, lb0, lb0, lb0, lb1, lb1,
		lb1, lb2, lb2, lb2, lb0, lb0, lb0, lb1, lb1, lb1, lb2, lb2, lb2 };

constexpr std::array<integer, NNEIGHBOR> zlb = { lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1,
		lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2, lb0, lb1, lb2 };

constexpr std::array<integer, NNEIGHBOR> xub = { ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub0, ub1, ub1, ub1, ub1, ub1,
		ub1, ub1, ub1, ub1, ub2, ub2, ub2, ub2, ub2, ub2, ub2, ub2, ub2 };

constexpr std::array<integer, NNEIGHBOR> yub = { ub0, ub0, ub0, ub1, ub1, ub1, ub2, ub2, ub2, ub0, ub0, ub0, ub1, ub1,
		ub1, ub2, ub2, ub2, ub0, ub0, ub0, ub1, ub1, ub1, ub2, ub2, ub2 };

constexpr std::array<integer, NNEIGHBOR> zub = { ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1,
		ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2, ub0, ub1, ub2 };

constexpr std::array<integer, NNEIGHBOR> bnd_size = { BW * BW * BW, BW * BW * NX, BW * BW * BW, BW * NX * BW, BW * NX
		* NX, BW * NX * BW, BW * BW * BW, BW * BW * NX, BW * BW * BW, NX * BW * BW, NX * BW * NX, NX * BW * BW, NX * NX
		* BW, NX * NX * NX, NX * NX * BW, NX * BW * BW, NX * BW * NX, NX * BW * BW, BW * BW * BW, BW * BW * NX, BW * BW
		* BW, BW * NX * BW, BW * NX * NX, BW * NX * BW, BW * BW * BW, BW * BW * NX, BW * BW * BW };

constexpr integer octlb0 = 0;
constexpr integer octlb1 = NX / 2 - 1;
constexpr integer octub0 = NX / 2;
constexpr integer octub1 = NX - 1;

constexpr std::array<integer, NCHILD> oct_xlb = { octlb0, octlb0, octlb0, octlb0, octlb1, octlb1, octlb1, octlb1 };

constexpr std::array<integer, NCHILD> oct_ylb = { octlb0, octlb0, octlb1, octlb1, octlb0, octlb0, octlb1, octlb1 };

constexpr std::array<integer, NCHILD> oct_zlb = { octlb0, octlb1, octlb0, octlb1, octlb0, octlb1, octlb0, octlb1 };

constexpr std::array<integer, NCHILD> oct_xub = { octub0, octub0, octub0, octub0, octub1, octub1, octub1, octub1 };

constexpr std::array<integer, NCHILD> oct_yub = { octub0, octub0, octub1, octub1, octub0, octub0, octub1, octub1 };

constexpr std::array<integer, NCHILD> oct_zub = { octub0, octub1, octub0, octub1, octub0, octub1, octub0, octub1 };

constexpr std::array<integer, NNEIGHBOR> dir_reverse = { 26, 25, 24, 23, 22, 21, 20, 19, 18,

17, 16, 15, 14, 13, 12, 11, 10, 9,

8, 7, 6, 5, 4, 3, 2, 1, 0 };

constexpr std::array<integer, NNEIGHBOR> dir_x = { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
		1, 1, 1, 1, 1, 1, 1 };

constexpr std::array<integer, NNEIGHBOR> dir_y = { -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1,
		-1, 0, 0, 0, 1, 1, 1 };

constexpr std::array<integer, NNEIGHBOR> dir_z = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1,
		-1, 0, 1, -1, 0, 1 };

integer ind3d(integer j, integer k, integer l) {
	return l + NX * (k + NX * j);
}

integer ind4d(integer p, integer j, integer k, integer l) {
	return l + NX * (k + NX * (j + NX * p));
}

std::size_t location_to_key(integer level, std::array<integer, NDIM> loc) {
	std::size_t key = 1;
	while (level > 0) {
		for (integer d = 0; d != NDIM; ++d) {
			key <<= 1;
			key |= (loc[d] & 1);
			loc[d] >>= 1;
		}
		level--;
	}
	return key;
}

bool location_is_phys_bnd(integer level, const std::array<integer, NDIM>& loc) {
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (dir_x[d] == +1 && loc[2] == (1 << level) - 1) {
			return true;
		}
		if (dir_x[d] == -1 && loc[2] == 0) {
			return true;
		}
		if (dir_y[d] == +1 && loc[1] == (1 << level) - 1) {
			return true;
		}
		if (dir_y[d] == -1 && loc[1] == 0) {
			return true;
		}
		if (dir_z[d] == +1 && loc[0] == (1 << level) - 1) {
			return true;
		}
		if (dir_z[d] == -1 && loc[0] == 0) {
			return true;
		}
	}
	return false;
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
		neighbor_status[d] = location_is_phys_bnd(level, location) ? COMPLETE : WAITING;
	}
}

hpx::future<std::array<real, PP * N3 / NCHILD>> node_server::get_expansions(integer c) const {
	return hpx::async([=]() {
		auto lptr = std::make_shared<std::array<real, PP * N3 / NCHILD>>();
		for (integer p = 0; p != PP; ++p) {
			for (integer j = oct_xlb[c]; j <= oct_xub[c]; ++j) {
				for (integer k = oct_ylb[c]; k <= oct_yub[c]; ++k) {
					for (integer l = oct_zlb[c]; l <= oct_zub[c]; ++l) {
						(*lptr)[ind4d(p, j - oct_xlb[c], k - oct_ylb[c], l - oct_zlb[c])] =
						L[ind4d(p,j,k,l)];
					}
				}
			}
		}
		return std::move(*lptr);
	});
}

hpx::future<std::array<real, PP * N3 / NCHILD>> node_server::get_multipoles() const {
	return hpx::async(hpx::launch::deferred, [=]() {
		auto mptr = std::make_shared<std::array<real, PP * N3 / NCHILD>>();
		real mc;
		for (integer p = 0; p != PP; ++p) {
			for (integer j = 0; j != NX; j += 2) {
				for (integer k = 0; k != NX; k += 2) {
					for (integer l = 0; l!= NX; l += 2) {
						mc = M[ind4d(p, j, k, l)] + M[ind4d(p, j, k, l + 1)]
						+ M[ind4d(p, j, k + 1, l)] + M[ind4d(p, j + 1, k, l)]
						+ M[ind4d(p, j, k + 1, l + 1)]
						+ M[ind4d(p, j + 1, k, l + 1)]
						+ M[ind4d(p, j + 1, k + 1, l)]
						+ M[ind4d(p, j + 1, k + 1, l + 1)];
						(*mptr)[ind4d(p, j >> 1, k >> 1, l >> 1)] = mc;
					}
				}
			}
		}
		return std::move(*mptr);
	});
}

void node_server::refine() {
	std::vector < hpx::future < hpx::id_type >> cids(NCHILD);
	for (integer ci = 0; ci != NCHILD; ++ci) {
		std::array<integer, NDIM> cloc(location);
		for (integer d = 0; d != NDIM; ++d) {
			cloc[d] <<= 1;
			cloc[d] += (ci >> d) & integer(1);
		}
		cids[ci] = hpx::new_ < node_server > (hpx::find_here(), my_id, level + 1, std::move(cloc));
	}
	is_leaf = false;
	hpx::wait_all (cids);
	for (integer ci = 0; ci != NCHILD; ++ci) {
		child_id[ci] = cids[ci].get();
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
		return bnd;
	});
}

void node_server::set_boundary(hpx::future<std::vector<real>> f, integer d) {
	neighbor_futures[d] = std::move(f);
	neighbor_status[d] = READY;
	input_condition.notify_one();
}

void node_server::set_multipoles(hpx::future<std::array<real, PP * N3 / NCHILD>> f, integer ci) {
	child_futures[ci] = std::move(f);
	child_status[ci] = READY;
	input_condition.notify_one();
}

void node_server::set_expansions(hpx::future<std::array<real, PP * N3 / NCHILD>> f) {
	parent_future = std::move(f);
	parent_status = READY;
	input_condition.notify_one();
}

void node_server::wait_for_signal() {
	hpx::this_thread::yield();

//	std::unique_lock<decltype(input_lock)> lock(input_lock);
//	input_condition.wait(lock);
}

void node_server::M2M(std::shared_ptr<std::array<real, PP * N3 / NCHILD>> mptr, integer c) {
	auto i = mptr->begin();
	for (integer p = 0; p != PP; ++p) {
		for (integer j = oct_xlb[c]; j <= oct_xub[c]; ++j) {
			for (integer k = oct_ylb[c]; k <= oct_yub[c]; ++k) {
				for (integer l = oct_zlb[c]; l <= oct_zub[c]; ++l) {
					M[ind4d(p, j, k, k)] = *i;
					L[ind4d(p, j, k, l)] = 0.0;
					++i;
				}
			}
		}
	}
}

void node_server::L2L(std::shared_ptr<std::array<real, PP * N3 / NCHILD>> l_coarse) {
	auto i = l_coarse->begin();
	for (integer p = 0; p != PP; ++p) {
		for (integer j = 0; j != NX; j += 2) {
			for (integer k = 0; k != NX; k += 2) {
				for (integer l = 0; l != NX; l += 2) {
					auto tmp = *i;
					L[ind4d(p, j, k, l)] += tmp;
					L[ind4d(p, j + 1, k, l)] += tmp;
					L[ind4d(p, j, k + 1, l)] += tmp;
					L[ind4d(p, j, k, l + 1)] += tmp;
					L[ind4d(p, j + 1, k + 1, l)] += tmp;
					L[ind4d(p, j + 1, k, l + 1)] += tmp;
					L[ind4d(p, j, k + 1, l + 1)] += tmp;
					L[ind4d(p, j + 1, k + 1, l + 1)] += tmp;
				}
			}
		}
	}
}

template<class Container>
void node_server::M2L(const Container& m, integer d) {
	const real delta = 1.0 / real(std::pow(integer(2), level));
	auto i = m.begin();
	for (integer j1 = xlb[d] + x_off[d]; j1 <= xub[d] + x_off[d]; ++j1) {
		const integer jlb = std::max(((j1 >> 1) - 1) << 1, integer(0));
		const integer jub = std::min(((j1 >> 1) + 1) << 1, NX - 1);
		for (integer k1 = ylb[d] + y_off[d]; k1 <= yub[d] + y_off[d]; ++k1) {
			const integer klb = std::max(((k1 >> 1) - 1) << 1, integer(0));
			const integer kub = std::min(((k1 >> 1) + 1) << 1, NX - 1);
			for (integer l1 = zlb[d] + z_off[d]; l1 <= zub[d] + z_off[d]; ++l1) {
				const integer llb = std::max(((l1 >> 1) - 1) << 1, integer(0));
				const integer lub = std::min(((l1 >> 1) + 1) << 1, NX - 1);
				for (integer p = 0; p != PP; ++p) {
					for (integer j0 = jlb; j0 <= jub; ++j0) {
						if (!is_leaf && std::abs(j0 - j1) < 2) {
							continue;
						}
						for (integer k0 = klb; k0 <= kub; ++k0) {
							if (!is_leaf && std::abs(k0 - k1) < 2) {
								continue;
							}
							for (integer l0 = llb; l0 <= lub; ++l0) {
								if (!is_leaf && std::abs(l0 - l1) < 2) {
									continue;
								}
								integer r2 = std::pow(j1 - j0, integer(2)) + std::pow(k1 - k0, integer(2))
										+ std::pow(l1 - l0, integer(2));
								if (r2 != 0) {
									const real rinv = 1.0 / (std::sqrt(real(r2)) * delta);
									L[ind4d(p, j0, k0, l0)] += *i * rinv;
								}
							}
						}
					}
				}
			}
		}
	}
}

std::vector<node_client> node_server::get_children_at_direction(integer d) const {
	std::vector<node_client> cids(cbnd_size[d]);
	auto i = cids.begin();
	for (integer j = dir_x[d] == +1 ? 1 : 0; j != (dir_x[d] == -1 ? 1 : 2); ++j) {
		for (integer k = dir_y[d] == +1 ? 1 : 0; k != (dir_y[d] == -1 ? 1 : 2); ++k) {
			for (integer l = dir_z[d] == +1 ? 1 : 0; l != (dir_z[d] == -1 ? 1 : 2); ++l) {
				*i = child_id[4 * j + k * 2 + l];
			}
		}
	}
	return cids;
}

void node_server::get_tree() {
	std::vector < hpx::future < hpx::id_type >> id_futs(NNEIGHBOR);
	std::vector < std::size_t > ids(NNEIGHBOR);
	auto i = ids.begin();
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (!location_is_phys_bnd(level, location)) {
			auto this_loc = location;
			this_loc[2] += dir_x[d];
			this_loc[1] += dir_y[d];
			this_loc[0] += dir_z[d];
			*i = location_to_key(level, std::move(this_loc));
			++i;
		}
	}
	auto f = hpx::find_ids_from_basename("fmmx_node", std::vector < std::size_t > (ids.begin(), i));
	hpx::wait_all(std::move(f));
	auto j = f.begin();
	for (integer d = 0; d != NNEIGHBOR; ++d) {
		if (!location_is_phys_bnd(level, location)) {
			neighbor_id[d] = j->get();
			++j;
		}
	}
}

void node_server::execute() {
	if (level < MAXLEVEL) {
		refine();
	}
	get_tree();

	printf("Executing thread at grid %i - %i %i %i\n", level, location[2], location[1], location[0]);
	while (1) {
		hpx::this_thread::sleep_for(boost::posix_time::milliseconds(1));
		if (false) {
			bool done;
			auto mul_fut = get_multipoles();
			integer j = 0;
			for (integer i = 0; i != NDIM; ++i) {
				j |= location[i] << i;
			}
			parent_id.set_multipoles(mul_fut, j);
			do {
				wait_for_signal();
				done = true;
				for (integer c = 0; c != NCHILD; ++c) {
					switch (child_status[c]) {
					case WAITING:
						done = false;
						break;
					case READY:
						M2M(std::make_shared<std::array<real, PP * N3 / NCHILD>>(child_futures[c].get()), c);
						child_status[c] = COMPLETE;
						break;
					default:
						break;
					}
				}
			} while (!done);
			for (integer d = 0; d != NNEIGHBOR; ++d) {
				if (d == center_dir) {
					continue;
				}
				auto bnd_fut = get_boundary(d);
				neighbor_id[d].set_boundary(bnd_fut, dir_reverse[d]);
			}
			M2L(M, center_dir);
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
						break;
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
					L2L(std::make_shared<std::array<real, PP * N3 / NCHILD>>(parent_future.get()));
					parent_status = COMPLETE;
					break;
				default:
					break;
				}
			} while (!done);
			if (!is_leaf) {
				for (integer ci = 0; ci != NCHILD; ++ci) {
					auto f = get_expansions(ci);
					child_id[ci].set_expansions(f);
				}
			}
		}
	}
}

node_server::node_server() :
		my_id(neighbor_id[center_dir]) {
	assert(false);
}

node_server::node_server(component_type* ptr) :
		base_type(ptr), my_id(neighbor_id[center_dir]) {
	node_client pid = hpx::naming::invalid_id;
	integer lev = 0;
	std::array<integer, NDIM> loc { { 0, 0, 0 } };
	initialize(std::move(pid), std::move(lev), std::move(loc));
}

node_server::node_server(component_type* ptr, node_client pid, integer lev, std::array<integer, NDIM> loc) :
		base_type(ptr), my_id(neighbor_id[center_dir]) {
	initialize(std::move(pid), std::move(lev), std::move(loc));
}

void node_server::initialize(node_client pid, integer lev, std::array<integer, NDIM> loc) {
	neighbors_set = false;
	location = loc;
	level = lev;
	parent_id = node_client(pid);
	my_id = get_gid();
	is_leaf = true;
	reset_children();
	reset_parent();
	reset_neighbors();
	key = location_to_key(level, std::move(loc));
	auto f = hpx::register_id_with_basename("fmmx_node", my_id, key);
	if (!f.get()) {
		printf("Failed to register id_with_basename\n");
		abort();
	}
	my_thread = hpx::thread([=]() {
		execute();
	});
}

node_server::~node_server() {
	my_thread.join();
}

