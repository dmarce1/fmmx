/*
 * main.cpp
 *
 *  Created on: Nov 9, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "node_server.hpp"

int hpx_main() {
	auto root = hpx::new_ < node_server > (hpx::find_here());
	return 0;
}
