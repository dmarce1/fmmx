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
#include <hpx/hpx_init.hpp>

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<node_server>, node_server);
typedef typename node_server::get_node_count_action get_node_count_action_t;
HPX_REGISTER_ACTION (get_node_count_action_t); //
typedef typename node_server::get_data_action get_data_action_t;
HPX_REGISTER_ACTION (get_data_action_t); //
typedef typename node_server::get_leaf_list_action get_leaf_list_action_t;
HPX_REGISTER_ACTION (get_leaf_list_action_t); //
typedef typename node_server::set_boundary_action set_boundary_action_t;
HPX_REGISTER_ACTION (set_boundary_action_t); //
typedef typename node_server::set_multipole_action set_multipole_action_t;
HPX_REGISTER_ACTION (set_multipole_action_t); //
typedef typename node_server::set_expansions_action set_expansions_action_t;
HPX_REGISTER_ACTION (set_expansions_action_t); //

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<silo_output>, silo_output);
typedef typename silo_output::do_output_action do_output_t;
HPX_REGISTER_ACTION (do_output_t);

int hpx_main() {
	auto silo = hpx::new_ < silo_output > (hpx::find_here());
	auto root = hpx::new_ < node_server > (hpx::find_here(), silo.get());
	auto f = hpx::register_id_with_basename("fmmx_node", root.get(), location_to_key(0, { 0, 0, 0 }));
	if (!f.get()) {
		printf("Failed to register id_with_basename\n");
		abort();
	}

	while (1) {
		hpx::this_thread::sleep_for(boost::posix_time::milliseconds(1000));
	}
	hpx::finalize();
	return 0;
}
