/*
 * node_client.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#include "node_client.hpp"
#include "node_server.hpp"

hpx::future<void> node_client::get_tree(std::vector<node_client> my_neighbors) {
	return hpx::async<typename node_server::get_tree_action>(id, std::move(my_neighbors));
}


bool node_client::operator==(const hpx::id_type& other) const {
	return id == other;
}
bool node_client::operator!=(const hpx::id_type& other) const {
	return id != other;
}

hpx::future<std::vector<node_client>> node_client::get_children() const {
	if (id != hpx::invalid_id) {
		return hpx::async<typename node_server::get_children_action>(id);
	} else {
		return hpx::make_ready_future(std::vector<node_client>(NCHILD, hpx::invalid_id));
	}
}

hpx::future<real> node_client::execute(real dt, integer rk) {
	return hpx::async<typename node_server::execute_action>(id, dt, rk);
}

hpx::future<void> node_client::refine(hpx::id_type myid) {
	return hpx::async<typename node_server::refine_action>(id, myid);
}

hpx::future<void> node_client::derefine() {
	assert(id != hpx::invalid_id);
	return hpx::async<typename node_server::derefine_action>(id, false);
}

hpx::future<void> node_client::destroy() {
	assert(id != hpx::invalid_id);
	return hpx::async<typename node_server::derefine_action>(id, true);
}

node_client::operator hpx::id_type() const {
	return id;
}

node_client::node_client() {
	id = hpx::naming::invalid_id;
}

node_client::node_client(const hpx::id_type& i) {
	id = i;
}

hpx::future<integer> node_client::get_node_count() const {
	return hpx::async<typename node_server::get_node_count_action>(id);

}

hpx::future<std::list<std::size_t>> node_client::get_leaf_list() const {
	return hpx::async<typename node_server::get_leaf_list_action>(id);
}

hpx::future<std::vector<double>> node_client::get_data() const {
	return hpx::async<typename node_server::get_data_action>(id);
}

node_client& node_client::operator=(const hpx::id_type& i) {
	id = i;
	return *this;
}

void node_client::set_multipoles(hpx::future<std::vector<real>>& f, integer ci) {
	if (id != hpx::naming::invalid_id) {
		hpx::apply<typename node_server::set_multipole_action>(id, std::move(f), ci);
	}
}

void node_client::set_expansions(hpx::future<std::vector<real>>& f) {
	if (id != hpx::invalid_id) {
		hpx::apply<typename node_server::set_expansions_action>(id, std::move(f));
	}
}

void node_client::set_boundary(hpx::future<std::vector<real>>& f, integer d) {
	if (id != hpx::invalid_id) {
		hpx::apply<typename node_server::set_boundary_action>(id, std::move(f), d);
	}
}

