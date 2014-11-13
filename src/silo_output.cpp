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
#include <thread>


template<class R, class ...Args>
R exec_on_separate_thread(R (*fptr)(Args...), Args ...args) {
	R data;
	std::thread([&](Args*...args) {
		data = (*fptr)(*args...);
	}, &args...).join();
	return data;

}


silo_output::silo_output() {
	//	printf("Silo in\n");
	received.resize((hpx::find_all_localities()).size());
	reset();
	//	printf("Silo out\n");

}
silo_output::~silo_output() {
}

void silo_output::do_output() {
	mutex0.lock();
#ifndef NO_OUTPUT
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
	exec_on_separate_thread(&DBPutZonelist2, db, "zones", nzones, int(NDIM), zone_nodes.data(), Nchild * nzones, 0, 0, 0,
			shapetype, shapesize, shapecnt, nshapes, olist);
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
	mutex0.unlock();
	reset();
}

void silo_output::reset() {
	mutex0.lock();
	current_index = 0;
	std::fill(received.begin(), received.end(), false);
	mutex0.unlock();
}

void silo_output::send_zones_to_silo(int proc_num_from, std::vector<zone> zones) {
	const int vertex_order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
	const int sz = zones.size();
	assert(!received[proc_num_from]);
	for (int i = 0; i < sz; i++) {
		silo_zone s;
		int j;
		s.fields = std::move(zones[i].fields);
		//	printf("%e %e %e %e %e %e\n", zones[i].position[0], zones[i].position[1], zones[i].position[2], zones[i].span[0], zones[i].span[1],
		//			zones[i].span[2]);
		for (int ci0 = 0; ci0 < Nchild; ci0++) {
			vertex v;
			int ci = vertex_order[ci0];
			for (int k = 0; k < NDIM; k++) {
				const double factor = (0.5 * double(2 * ((ci >> k) & 1) - 1));
				v[k] = zones[i].position[k] + zones[i].span[k] * factor;
			}
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
	mutex0.lock();
	received[proc_num_from] = true;
	printf("Receive output from %i\n", proc_num_from);
	if (std::all_of(received.begin(), received.end(), [](bool b) {return b;})) {
		printf("Doing output\n");
		mutex0.unlock();
		do_output();
		mutex0.lock();
	}
	mutex0.unlock();

}
