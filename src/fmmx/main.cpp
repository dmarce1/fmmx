/*
 * main.cpp
 *
 *  Created on: Nov 9, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "silo_output.hpp"
#include "node_server.hpp"
#include "key.hpp"
#include <chrono>

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<node_server>, node_server);
typedef typename node_server::get_node_count_action get_node_count_action_t;
HPX_REGISTER_ACTION(get_node_count_action_t);
//
typedef typename node_server::get_tree_action get_tree_action_t;
HPX_REGISTER_ACTION(get_tree_action_t);

typedef typename node_server::get_children_action get_children_action_t;
HPX_REGISTER_ACTION(get_children_action_t);
//
typedef typename node_server::get_data_action get_data_action_t;
HPX_REGISTER_ACTION(get_data_action_t);
//
typedef typename node_server::get_leaf_list_action get_leaf_list_action_t;
HPX_REGISTER_ACTION(get_leaf_list_action_t);
//
typedef typename node_server::set_boundary_action set_boundary_action_t;
HPX_REGISTER_ACTION(set_boundary_action_t);
//
typedef typename node_server::set_multipole_action set_multipole_action_t;
HPX_REGISTER_ACTION(set_multipole_action_t);
//
typedef typename node_server::set_expansions_action set_expansions_action_t;
HPX_REGISTER_ACTION(set_expansions_action_t);
//
typedef typename node_server::execute_action execute_action_t;
HPX_REGISTER_ACTION(execute_action_t);
//

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<silo_output>, silo_output);
typedef typename silo_output::do_output_action do_output_t;
HPX_REGISTER_ACTION(do_output_t);

int hpx_main() {
	std::chrono::time_point<std::chrono::system_clock> start, end;

	node_client root_client;

	auto silo = hpx::new_<silo_output>(hpx::find_here());

	auto sout = silo.get();

	auto root = hpx::new_<node_server>(hpx::find_here(), sout);

	root_client = root.get();

	auto f = hpx::register_id_with_basename("fmmx_node", root_client, location_to_key(0, std::array<integer, NDIM> { {
			0, 0, 0 } }));

	if (!f.get()) {
		printf("Failed to register id_with_basename\n");
		abort();
	}

	start = std::chrono::system_clock::now();
	root_client.refine(root_client).get();
	root_client.get_tree().get();
	end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "finished refinement in " << elapsed_seconds.count() << "s\n";

	//++fnum;

	real tmax = 0.1;
	start = std::chrono::system_clock::now();
	integer fnum = 0;
	printf( "Executing...\n");
	real dt = root_client.execute(0.0,0).get();
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	std::cout << "finished computation in " << elapsed_seconds.count() << "s\n";
	root_client.execute(dt,1).get();
	real t = 0.0;
	for( integer z = 0; z != 25; z++) {
	//while (t < tmax) {
		t += dt;
		printf("%e %e\n", double(t), double(dt));
		dt = root_client.execute(dt,0).get();
				root_client.execute(dt,1).get();
	}
	root_client.execute(dt,0).get();

	std::list<std::size_t> leaf_list = root_client.get_leaf_list().get();
	printf("%li leaves detected by root\n", leaf_list.size());
	auto f1 = hpx::async<typename silo_output::do_output_action>(sout, std::move(leaf_list), 0);
	f1.get();


	auto f0 = root_client.destroy();
	f0.get();
	hpx::finalize();
	return 0;
}