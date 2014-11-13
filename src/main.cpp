/*
 * main.cpp
 *
 *  Created on: Nov 9, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "node_server.hpp"
#include <hpx/hpx_init.hpp>


HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<node_server>, node_server);

int hpx_main() {
	auto root = hpx::new_ < node_server > (hpx::find_here());
	while(1) {
		hpx::this_thread::sleep_for(boost::posix_time::milliseconds(1000));
	}
	hpx::finalize();
	return 0;
}
