/*
 * node_server.hpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#ifndef __NODE_SERVER__HPP
#define __NODE_SERVER__HPP

#include <hpx/lcos/local/condition_variable.hpp>

class node_server: public hpx::components::managed_component_base<node_server, hpx::components::detail::this_type,
		hpx::traits::construct_with_back_ptr> {
public:
	using component_type = hpx::components::managed_component<node_server>;
	using base_type = hpx::components::managed_component_base<node_server, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr>;
private:
	std::array<real, PP * N3> M;
	std::array<real, PP * N3> L;
	node_client parent_id;
	std::array<node_client, NCHILD> child_id;
	std::array<node_client, NNEIGHBOR> neighbor_id;
	std::atomic<integer> parent_status;
	std::array<std::atomic<integer>, NCHILD> child_status;
	std::array<std::atomic<integer>, NNEIGHBOR> neighbor_status;
	std::array<hpx::future<std::vector<real>>, NNEIGHBOR> neighbor_futures;
	std::array<hpx::future<std::array<real, PP * N3 / NCHILD>>, NCHILD> child_futures;
	hpx::future<std::array<real, PP * N3 / NCHILD>> parent_future;
	hpx::lcos::local::condition_variable input_condition;
	hpx::lcos::local::spinlock input_lock;
	std::array<integer, NDIM> location;
	integer level;
	hpx::thread my_thread;
	bool is_leaf;
	void initialize(node_client, integer, std::array<integer, NDIM>);
	node_client& my_id;
public:
	node_server(component_type*);
	node_server(component_type*, node_client, integer, std::array<integer, NDIM>);
	~node_server();
	void get_tree();
	hpx::future<std::array<real, PP * N3 / NCHILD>> get_multipoles() const;
	hpx::future<std::array<real, PP * N3 / NCHILD>> get_expansions(integer ci) const;
	hpx::future<std::vector<real>> get_boundary(integer d) const;
	void set_boundary(hpx::future<std::vector<real>> f, integer d);
	void set_multipoles(hpx::future<std::array<real, PP * N3 / NCHILD>> f, integer ci);
	void set_expansions(hpx::future<std::array<real, PP * N3 / NCHILD>>);
	std::vector<node_client> get_children_at_direction(integer) const;
	void wait_for_signal();
	void M2M(std::shared_ptr<std::array<real, PP * N3 / NCHILD>>, integer);
	template<class Container>
	void M2L(const Container&, integer);
	void L2L(std::shared_ptr<std::array<real, PP * N3 / NCHILD>>);
	void reset_children();
	void reset_parent();
	void reset_neighbors();
	void execute();
	void refine();
	//
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_boundary, set_boundary_action);
	//
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_multipoles, set_multipole_action);
	//
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_expansions, set_expansions_action);
	//
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_children_at_direction, get_children_at_direction_action);
};

#endif

