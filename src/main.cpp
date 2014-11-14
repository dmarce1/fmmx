/*
 * main.cpp
 *
 *  Created on: Nov 9, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "silo_output.hpp"
#include "node_server.hpp"
#include <hpx/hpx_init.hpp>

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<node_server>, node_server);
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<silo_output>, silo_output);

int hpx_main() {
	auto silo = hpx::new_ < silo_output > (hpx::find_here());
	auto root = hpx::new_ < node_server > (hpx::find_here(), silo.get());
	while (1) {
		hpx::this_thread::sleep_for(boost::posix_time::milliseconds(1000));
	}
	hpx::finalize();
	return 0;
}
