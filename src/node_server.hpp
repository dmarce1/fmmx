/*
 * node_server.hpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#ifndef __NODE_SERVER__HPP
#define __NODE_SERVER__HPP

#include <hpx/lcos/local/counting_semaphore.hpp>
#include <boost/serialization/list.hpp>

class node_server: public hpx::components::managed_component_base<node_server> {
private:
	static hpx::id_type output;
	std::vector<real> M;
	std::vector<real> L;
	node_client parent_id;
	std::array<node_client, NCHILD> child_id;
	std::array<node_client, NNEIGHBOR> neighbor_id;
	std::atomic<integer> parent_status;
	std::atomic<bool> neighbors_set;
	std::array<std::atomic<integer>, NCHILD> child_status;
	std::array<std::atomic<integer>, NNEIGHBOR> neighbor_status;
	std::array<hpx::future<std::vector<real>>, NNEIGHBOR> neighbor_futures;
	std::array<hpx::future<std::vector<real>>, NCHILD> child_futures;
	hpx::future<std::vector<real>> parent_future;
	mutable hpx::lcos::local::counting_semaphore input_condition;
	mutable hpx::lcos::local::counting_semaphore data_ready;
	std::array<integer, NDIM> location;
	hpx::id_type my_id;
	real dx;
	integer level;
	hpx::thread my_thread;
	bool is_leaf;
	std::uint64_t key;
	void initialize(node_client, integer, std::array<integer, NDIM>);
public:
	node_server();
	node_server(hpx::id_type);
	node_server(node_client, integer, std::array<integer, NDIM>);
	~node_server();
	void get_tree();
	void init_t0();
	integer get_node_count() const;
	std::vector<real> get_data() const;
	std::list<std::size_t> get_leaf_list() const;
	hpx::future<std::vector<real>> get_multipoles() const;
	hpx::future<std::vector<real>> get_expansions(integer ci) const;
	hpx::future<std::vector<real>> get_boundary(integer d) const;
	void set_boundary(hpx::future<std::vector<real>> f, integer d);
	void set_multipoles(hpx::future<std::vector<real>> f, integer ci);
	void set_expansions(hpx::future<std::vector<real>>);
	void wait_for_signal() const;
	void M2M(const std::vector<real>&, integer);
	void M2L(const std::vector<real>&, integer);
	void L2L(const std::vector<real>&);
	void reset_children();
	void reset_parent();
	void reset_neighbors();
	void execute();
	void refine();
	//
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_node_count, get_node_count_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_data, get_data_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_leaf_list, get_leaf_list_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_boundary, set_boundary_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_multipoles, set_multipole_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_expansions, set_expansions_action);
	//
	//
};

#endif

