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

#include "math.hpp"
#include "silo_output.hpp"
#include "node_client.hpp"
#include "key.hpp"
#include <thread>

silo_output::silo_output() {
}
/*
 template<class R, class ...Args1, class ...Args2>
 R exec_on_separate_thread(R (*fptr)(Args1...), Args2 ...args) {
 R data;
 std::thread([&](Args2*...args) {
 data = (*fptr)(*args...);
 }, &args...).join();
 return data;

 }
 */
void silo_output::do_output(std::list<std::size_t> node_list, integer filenum) {
	current_index = 0;
#ifndef NO_OUTPUT
	std::vector<std::size_t> id_list(node_list.begin(), node_list.end());
	auto id_list_ptr = &id_list;
	std::vector<hpx::future<hpx::id_type> > node_futs = hpx::find_ids_from_basename("fmmx_node", id_list);
	std::vector<hpx::future<void>> data_futs(node_futs.size());
	for (integer i0 = 0; i0 != node_futs.size(); ++i0) {
		data_futs[i0] = (node_client(node_futs[i0].get()).get_data()).then(hpx::util::unwrapped([=](std::vector<double> data) {
			auto iter = data.begin();
			auto key = (*id_list_ptr)[i0];
			constexpr integer vertex_order[8] = {0, 1, 3, 2, 4, 5, 7, 6};
			std::array<integer, NDIM> loc;
			std::array<double, NDIM> corner;
			integer lev;
			double span;
			key_to_location(key, &lev, &loc);
			for( integer a = 0; a != NDIM; ++a) {
				corner[a] = double(loc[a]) / double(std::pow(2,lev));
			}
			integer cnt = 0;
			span = 1.0 / double(NX*std::pow(2,lev));
			for (integer j0 = 0; j0 != SILO_NX; ++j0) {
				for (integer k0 = 0; k0 != SILO_NX; ++k0) {
					for (integer l0 = 0; l0 != SILO_NX; ++l0) {
						const integer jp = j0 / HYDRO_P;
						const integer kp = k0 / HYDRO_P;
						const integer lp = l0/ HYDRO_P;
						const integer jc = j0 % HYDRO_P;
						const integer kc = k0 % HYDRO_P;
						const integer lc = l0 % HYDRO_P;
						silo_zone s;
						int j;
						for( integer k = 0; k != HYDRO_NF; ++k) {
							s.fields[k] = *iter;
							iter++;
							cnt++;
						}
						auto dx0 = span * gauss_weights[HYDRO_P-1][jc] / real(2);
						auto dy0 = span * gauss_weights[HYDRO_P-1][kc] / real(2);
						auto dz0 = span * gauss_weights[HYDRO_P-1][lc] / real(2);
						auto x0 = corner[0] + (span/real(2))*(gauss_vol_lb[HYDRO_P-1][jc]+real(1));
						auto y0 = corner[1] + (span/real(2))*(gauss_vol_lb[HYDRO_P-1][kc]+real(1));
						auto z0 = corner[2] + (span/real(2))*(gauss_vol_lb[HYDRO_P-1][lc]+real(1));
						for (int ci0 = 0; ci0 < Nchild; ci0++) {
							vertex v;
							int ci = vertex_order[ci0];
							v[0] = real(jp)*span + real((ci >> 0) & 1)*dx0 + x0;
							v[1] = real(kp)*span + real((ci >> 1) & 1)*dy0 + y0;
							v[2] = real(lp)*span + real((ci >> 2) & 1)*dz0 + z0;
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

	hpx::wait_all(data_futs);

	std::thread([=]() {
		constexpr int nshapes = 1;
		const int nnodes = nodedir.size();
		const int nzones = zonedir.size();
		int shapesize[1] = {Nchild};
		int shapetype[1] = {DB_ZONETYPE_HEX};
		int shapecnt[1] = {nzones};
		DBfile * db;
		DBoptlist * olist;
		std::vector<double> coord_vectors[NDIM];
		std::string coordname_strs[NDIM];
		double* coords[NDIM];
		const char* coordnames[NDIM];

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

	auto this_sz = zonedir.size();
	std::vector<int> zone_nodes(this_sz * Nchild);
	int zni = 0;
	for (auto zi = zonedir.begin(); zi != zonedir.end(); zi++) {
		for (int ci0 = 0; ci0 < Nchild; ci0++) {
			zone_nodes[zni] = zi->vertices[ci0];
			zni++;
		}
	}

	olist = DBMakeOptlist( 1);
	char* fname;
	asprintf(&fname, "X.%i.silo", filenum);
	db = DBCreateReal( fname, DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
	DBPutZonelist2( db, "zones", nzones, int(NDIM), zone_nodes.data(), Nchild * nzones, 0, 0,
			0, shapetype, shapesize, shapecnt, nshapes, olist);
	DBPutUcdmesh( db, "mesh", int(NDIM), const_cast<char**>(coordnames),
			reinterpret_cast<void*>(coords), nnodes, nzones, "zones", static_cast<const char*>(nullptr),
			(int) DB_DOUBLE, olist);

	std::vector<double> data(nzones);
	std::array<char[32], HYDRO_NF> field_names;
	sprintf(field_names[0], "D");
	sprintf(field_names[1], "SX");
	sprintf(field_names[2], "SY");
	sprintf(field_names[3], "SZ");
	sprintf(field_names[4], "E");
	for (int fi = 0; fi != HYDRO_NF; ++fi) {
		int i = 0;
		for (auto zi = zonedir.begin(); zi != zonedir.end(); zi++) {
			data[i] = zi->fields[fi];
			i++;
		}
		DBPutUcdvar1( db, const_cast<const char*>(field_names[fi]), "mesh",
				reinterpret_cast<void*>(data.data()), nzones, static_cast<void*>(nullptr), 0, (int) DB_DOUBLE,
				(int) DB_ZONECENT, olist);
	}

	DBClose(db);
	free(fname);
	zonedir.clear();
	nodedir.clear();
	printf("Output Done\n");
})	.join();
#endif
}
