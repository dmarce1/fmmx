/*
 * node_server.hpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#ifndef __NODE_SERVER__HPP
#define __NODE_SERVER__HPP

#include "defs.hpp"
#include "hydro.hpp"
#include <boost/serialization/list.hpp>
#include <boost/atomic.hpp>

class node_server: public hpx::components::managed_component_base<node_server> {
private:
	static hpx::id_type output;
	std::vector<real> M;
	std::vector<real> L;
	std::vector<real> phi;
	std::vector<real> gx, gy, gz;
	node_client parent_id;
	std::array<node_client, NCHILD> child_id;
	std::array<node_client, NNEIGHBOR> neighbor_id;
	std::array<bool,NNEIGHBOR> neighbor_is_leaf;
	mutable hpx::lcos::local::spinlock L_lock;
	mutable hpx::lcos::local::spinlock M_lock;
	std::array<integer, NDIM> location;
	hpx::id_type my_id;
	real dx;
	integer level;
	bool is_leaf;
	std::uint64_t key;
	void initialize(node_client, integer, std::array<integer, NDIM>);
	void reset();

	bool amr_bnd;

	boost::atomic<integer> child_done_cnt;
	boost::atomic<integer> neighbor_done_cnt;
	integer step_cnt;
	std::pair<hpx::promise<void>, hpx::promise<void>> exe_pair;
	hpx::promise<void>* exe_promise;

	hydro_vars hydro_state;
	real this_dt;
	bool hydro_updated;

public:
	bool refine_me() const;
	node_server();
	node_server(hpx::id_type);
	node_server(node_client, integer, std::array<integer, NDIM>);
	~node_server();
	void get_tree(
			std::vector<node_client> my_neighbors = std::vector<node_client>(NNEIGHBOR, node_client(hpx::invalid_id)));
	void init_t0();
	integer get_node_count() const;
	std::vector<double> get_data() const;
	std::list<std::size_t> get_leaf_list() const;
	void send_multipoles();
	void send_expansions();
	std::vector<node_client> get_children() const;
	hpx::future<std::vector<real>> get_multipoles();
	hpx::future<std::vector<real>> get_expansions(integer ci) const;
	hpx::future<std::vector<real>> get_boundary(integer d) const;
	void set_boundary(hpx::future<std::vector<real>> f, integer d);
	void set_multipoles(hpx::future<std::vector<real>> f, integer ci);
	void set_expansions(hpx::future<std::vector<real>>);
	void wait_for_signal() const;
	void M2M(const std::vector<real>&, integer);
	void M2L(const std::vector<real>&, integer, integer);
	void L2L(const std::vector<real>&);
	void reset_parent();
	real execute(real, integer);
	void derefine(bool);
	void refine(hpx::id_type id);
	//
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_tree, get_tree_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_children, get_children_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, derefine, derefine_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, refine, refine_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, execute, execute_action); //
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

