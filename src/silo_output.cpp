/*
 * silo_output.cpp
 *
 *  Created on: Nov 13, 2014
 *      Author: dmarce1
 */

/*
 * silo_output.hpp
 *
 *  Created on: May 31, 2014
 *      Author: dmarce1
 */

#include "silo_output.hpp"
#include "node_client.hpp"
#include "key.hpp"
#include <thread>

#include <hpx/lcos/when_any.hpp>

template<class R, class ...Args>
R exec_on_separate_thread(R (*fptr)(Args...), Args ...args) {
	R data;
	std::thread([&](Args*...args) {
		data = (*fptr)(*args...);
	}, &args...).join();
	return data;

}

void silo_output::do_output(std::list<std::size_t> node_list) {

	std::vector<hpx::future<hpx::id_type> >
	find_ids_from_basename(char const * base_name, std::vector<std::size_t> const & ids);

	std::vector<std::size_t> id_list(node_list.begin(), node_list.end());
	auto id_list_ptr = &id_list;
	std::vector<hpx::future<hpx::id_type> > node_futs = hpx::find_ids_from_basename("fmmx_node", std::move(id_list));
	std::vector<hpx::future<void>> data_futs(node_futs.size());

#ifndef NO_OUTPUT
	for (integer i0 = 0; i0 != node_futs.size(); ++i0) {
		data_futs[i0] = (node_client(node_futs[i0].get()).get_data()).then(
				hpx::util::unwrapped([=](std::vector<real> data) {
					auto iter = data.begin();
					auto key = (*id_list_ptr)[i0];
					constexpr integer vertex_order[8] = {0, 1, 3, 2, 4, 5, 7, 6};
					std::array<integer, NDIM> loc;
					std::array<real, NDIM> corner;
					integer lev;
					real span;
					key_to_location(key, &lev, &loc);
					for( integer a = 0; a != NDIM; ++a) {
						corner[a] = real(loc[a]) / real(std::pow(2,lev));
					}
					span = 1.0 / real(NX*std::pow(2,lev));
					for (integer j0 = 0; j0 != NX; ++j0) {
						for (integer k0 = 0; k0 != NX; ++k0) {
							for (integer l0 = 0; l0 != NX; ++l0) {
								silo_zone s;
								int j;
								for( integer k = 0; k != NF; ++k) {
									s.fields[k] = *iter;
									iter++;
								}
								for (int ci0 = 0; ci0 < Nchild; ci0++) {
									vertex v;
									int ci = vertex_order[ci0];
									v[0] = (j0 + (0.5 * real(2 * ((ci >> 0) & 1) - 1)))*span + corner[0];
									v[1] = (k0 + (0.5 * real(2 * ((ci >> 1) & 1) - 1)))*span + corner[1];
									v[2] = (l0 + (0.5 * real(2 * ((ci >> 2) & 1) - 1)))*span + corner[2];
									mutex0.lock();
									auto iter = nodedir.find(v);
									if (iter == nodedir.end()) {
										j = current_index;
										v.index = j;
										current_index++;
										nodedir.insert(std::move(v));
									} else {
										j = iter->index;
									}
									mutex0.unlock();
									s.vertices[ci0] = j;
								}
								mutex0.lock();
								zonedir.push_back(std::move(s));
								mutex0.unlock();
							}
						}
					}
				}));
	}

	hpx::wait_all(std::move(data_futs));
	constexpr int nshapes = 1;
	const int nnodes = nodedir.size();
	const int nzones = zonedir.size();
	int shapesize[1] = { Nchild };
	int shapetype[1] = { DB_ZONETYPE_HEX };
	int shapecnt[1] = { nzones };
	DBfile * db;
	DBoptlist * olist;
	std::vector<double> coord_vectors[NDIM];
	std::string coordname_strs[NDIM];
	double* coords[NDIM];
	const char* coordnames[NDIM];
	const int nfields = (zonedir.begin())->fields.size();

	for (int di = 0; di < NDIM; di++) {
		coord_vectors[di].resize(nnodes);
		coords[di] = coord_vectors[di].data();
		coordname_strs[di] = ('x' + char(di));
		coordnames[di] = coordname_strs[di].c_str();
	}

	for (auto ni = nodedir.begin(); ni != nodedir.end(); ni++) {
		for (int di = 0; di < NDIM; di++) {
			coords[di][ni->index] = (*ni)[di];
			//	printf("%e ", (*ni)[di]);
		}
		//	printf("\n");
	}
	nodedir.clear();

	std::vector<int> zone_nodes(zonedir.size() * Nchild);
	int zni = 0;
	for (auto zi = zonedir.begin(); zi != zonedir.end(); zi++) {
		for (int ci0 = 0; ci0 < Nchild; ci0++) {
			zone_nodes[zni] = zi->vertices[ci0];
			zni++;
		}
	}

	olist = exec_on_separate_thread(&DBMakeOptlist, 1);
	db = exec_on_separate_thread(&DBCreateReal, "X.silo", DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
	exec_on_separate_thread(&DBPutZonelist2, db, "zones", nzones, int(NDIM), zone_nodes.data(), Nchild * nzones, 0, 0,
			0, shapetype, shapesize, shapecnt, nshapes, olist);
	exec_on_separate_thread(&DBPutUcdmesh, db, "mesh", int(NDIM), const_cast<char**>(coordnames),
			reinterpret_cast<DB_DTPTR2>(coords), nnodes, nzones, "zones", static_cast<const char*>(nullptr),
			(int) DB_DOUBLE, olist);

	std::vector<double> data(nzones);
	char fname[2];
	fname[1] = '\0';
	for (int fi = 0; fi != nfields; ++fi) {
		int i = 0;
		for (auto zi = zonedir.begin(); zi != zonedir.end(); zi++) {
			data[i] = zi->fields[fi];
			i++;
		}
		fname[0] = 'A' + char(fi);
		exec_on_separate_thread(&DBPutUcdvar1, db, const_cast<const char*>(fname), "mesh",
				reinterpret_cast<DB_DTPTR1>(data.data()), nzones, static_cast<DB_DTPTR1>(nullptr), 0, (int) DB_DOUBLE,
				(int) DB_ZONECENT, olist);
	}

	zonedir.clear();

	exec_on_separate_thread(DBClose, db);
#endif
}
