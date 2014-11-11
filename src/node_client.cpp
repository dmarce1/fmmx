/*
 * node_client.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "node_server.hpp"

node_client::operator hpx::id_type() const {
	return id;
}

node_client::node_client() {
	id = hpx::naming::invalid_id;
}

node_client::node_client(const hpx::id_type& i) {
	id = i;
}

void node_client::set_neighbors(std::vector<node_client> ns) {
	if (id != hpx::naming::invalid_id) {
		hpx::apply<typename node_server::set_neighbors_action>(id, std::move(ns));
	}

}

hpx::future<std::vector<node_client>> node_client::get_children_at_direction(integer d) {
	if (id != hpx::naming::invalid_id) {
		return hpx::async<typename node_server::get_children_at_direction_action>(id, d);
	} else {
		std::vector<node_client> tmp;
		return hpx::make_ready_future(tmp);
	}
}

node_client& node_client::operator=(const hpx::id_type& i) {
	id = i;
	return *this;
}

void node_client::set_multipoles(hpx::future<std::array<real, PP * N3 / NCHILD>>& f, integer ci) {
	if (id != hpx::naming::invalid_id) {
		hpx::apply<typename node_server::set_multipole_action>(id, std::move(f), ci);
	}
}

void node_client::set_expansions(hpx::future<std::array<real, PP * N3 / NCHILD>>& f) {
	if (id != hpx::invalid_id) {
		hpx::apply<typename node_server::set_expansions_action>(id, std::move(f));
	}
}

void node_client::set_boundary(hpx::future<std::vector<real>>& f, integer d) {
	if (id != hpx::invalid_id) {
		hpx::apply<typename node_server::set_boundary_action>(id, std::move(f), d);
	}
}

