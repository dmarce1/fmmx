/*
 * main.cpp
 *
 *  Created on: Nov 9, 2014
 *      Author: dmarce1
 */
#include <mpi.h>
#include "node_client.hpp"
#include "silo_output.hpp"
#include "node_server.hpp"
#include "key.hpp"
#include <chrono>
#ifndef NDEBUG
#include <fenv.h>
#endif

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<node_server>, node_server);
typedef typename node_server::get_node_count_action get_node_count_action_t;
HPX_REGISTER_ACTION(get_node_count_action_t);

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<silo_output>, silo_output);
typedef typename silo_output::do_output_action do_output_t;
HPX_REGISTER_ACTION(do_output_t);

extern void hydro_initialize();

int hpx_main() {

	hydro_initialize();

#ifndef NDEBUG
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
#endif

	std::chrono::time_point<std::chrono::system_clock> start, end;

	node_client root_client;

	auto silo = hpx::new_<silo_output>(hpx::find_here());

	auto sout = silo.get();

	auto root = hpx::new_<node_server>(hpx::find_here(), sout);

	root_client = root.get();

	auto f = hpx::register_id_with_basename("fmmx_node", root_client, location_to_key(0, std::array<integer, NDIM> { { 0, 0, 0 } }));

	if (!f.get()) {
		printf("Failed to register id_with_basename\n");
		abort();
	}

	start = std::chrono::system_clock::now();

	root_client.set_me(root_client).get();
	for (integer l = 0; l < MAXLEVEL; ++l) {
		root_client.refine().get();
		root_client.get_tree().get();
		printf("Refined to level %li\n", l + 1);
	}

	//++fnum;
	std::list<std::size_t> leaf_list = root_client.get_leaf_list().get();
	printf("%li leaves detected by root\n", leaf_list.size());

	real tmax = 100.0;
	printf("Executing...\n");
	real dt;
	integer step = 0;
	real t = real(0);
	dt = real(0);
	root_client.hydro_project(0).get();
	root_client.hydro_restrict(0).get();
	root_client.hydro_amr_prolong(0).get();
	root_client.hydro_exchange(0,HYDRO_BND).get();
	root_client.hydro_exchange(0,GRAV_BND).get();
	root_client.execute(0).get();
	root_client.hydro_exchange(0,GRAV_BND).get();
	auto f1 = hpx::async<typename silo_output::do_output_action>(sout, leaf_list, 0);
	f1.get();
	const real dt_output = 1.0e-1;
	integer out_cnt = 1;
	while (t < tmax) {
	//	break;
		real tstart = MPI_Wtime();
		for (integer rk = 0; rk != HYDRO_RK; ++rk) {

			auto tfut = root_client.hydro_next_du(rk);
			if (rk == 0) {
				dt = tfut.get().first / real(2 * HYDRO_P + 1) * cfl[HYDRO_RK - 1];
				printf("%li %e %e\n", step, double(t), double(dt));
			} else {
				tfut.get();
			}

			root_client.hydro_next_u(rk, dt).get();
			const integer rk0 = (rk != HYDRO_RK - 1 ? rk + 1 : 0);
			root_client.hydro_amr_prolong(rk0).get();
			root_client.hydro_exchange(rk0,HYDRO_BND).get();
			root_client.hydro_project(rk0).get();
			root_client.hydro_restrict(rk0).get();
			root_client.hydro_amr_prolong(rk0).get();
			root_client.hydro_exchange(rk0,HYDRO_BND).get();
			root_client.execute(rk0).get();
			root_client.hydro_exchange(rk0,GRAV_BND).get();
		}
		++step;
		if (out_cnt*dt_output < t) {
			f1 = hpx::async<typename silo_output::do_output_action>(sout, leaf_list, out_cnt++);
			f1.get();
		}
		t += dt;
		real tend = MPI_Wtime();
//		printf( "%e\n", tend - tstart);
	}
//	 leaf_list = root_client.get_leaf_list().get();

	/*real t = 0.0;
	 real dt = 0.0;
	 for (integer z = 0; z != 50; z++) {
	 //	while (t < tmax) {
	 dt = root_client.execute(0).get();
	 for (integer rk = 1; rk < HYDRO_RK; ++rk) {
	 root_client.execute(rk).get();
	 }
	 t += dt;
	 leaf_list = root_client.get_leaf_list().get();
	 }*/
//root_client.execute(dt, 0).get();
//printf("%li leaves detected by root\n", leaf_list.size());
	auto f0 = root_client.destroy();
	f0.get();
	hpx::finalize();
	return 0;
}
