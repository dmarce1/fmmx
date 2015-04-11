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
#include "rk.hpp"
#include <boost/serialization/list.hpp>
#include <boost/atomic.hpp>



class node_server: public hpx::components::managed_component_base<node_server> {
private:
	static hpx::id_type output;
	integer this_rk;
	std::vector<real> M;
	std::vector<real> L;
//	std::vector<real> phi;
//	std::vector<real> gx, gy, gz;
	node_client parent_id;
	std::array<node_client, NCHILD> child_id;
	std::array<node_client, NNEIGHBOR> neighbor_id;
	std::array<bool, NNEIGHBOR> neighbor_is_leaf;
	mutable hpx::lcos::local::spinlock L_lock;
	mutable hpx::lcos::local::spinlock M_lock;
	mutable hpx::lcos::local::spinlock R_lock;
	std::array<integer, NDIM> location;
	hpx::id_type my_id;
	real dx;
	integer level;
	bool is_leaf;
	std::uint64_t key;
	void initialize(node_client, integer, std::array<integer, NDIM>);
	void reset();

	boost::atomic<integer> fmm_child_done_cnt;
	boost::atomic<integer> fmm_neighbor_done_cnt;
	integer step_cnt;
	std::pair<hpx::promise<void>, hpx::promise<void>> exe_pair;
	hpx::promise<void>* exe_promise;

	std::shared_ptr<hydro> hydro_vars;
public:

	bool refine_me() const;
	node_server();
	node_server(hpx::id_type);
	node_server(node_client, integer, std::array<integer, NDIM>);
	~node_server();
	void get_tree(std::vector<node_client> my_neighbors = std::vector<node_client>(NNEIGHBOR, node_client(hpx::invalid_id)));
	void init_t0();
	integer get_node_count() const;
	std::vector<double> get_data() const;
	std::list<std::size_t> get_leaf_list() const;
	void send_multipoles();
	void send_expansions();
	std::vector<node_client> get_children() const;
	hpx::future<std::vector<real>> get_multipoles();
	hpx::future<std::vector<real>> get_expansions(integer ci) const;
	hpx::future<std::vector<real>> get_fmm_boundary(integer d) const;
	void set_fmm_boundary(hpx::future<std::vector<real>> f, integer d);
	void set_multipoles(hpx::future<std::vector<real>> f, integer ci);
	void set_expansions(hpx::future<std::vector<real>>);
	void M2M(const std::vector<real>&, integer);
	void M2L(const std::vector<real>&, integer, integer);
	void L2L(const std::vector<real>&);
	real execute(integer rk);
	void derefine(bool);
	bool child_is_amr(integer ci, integer dir) const;
	bool is_amr(integer dir) const;
	bool refine();
	void refine_proper();
	void set_me(hpx::id_type id);
	bool is_phys_bound(integer dir) const;
	//
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_tree, get_tree_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_children, get_children_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, derefine, derefine_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, refine, refine_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, refine_proper, refine_proper_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_me, set_me_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, execute, execute_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_node_count, get_node_count_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_data, get_data_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_leaf_list, get_leaf_list_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_fmm_boundary, set_fmm_boundary_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_multipoles, set_multipole_action); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, set_expansions, set_expansions_action);
	//

	void hydro_exchange(integer, exchange_type); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, hydro_exchange, hydro_exchange_action);

	void hydro_amr_prolong(integer); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, hydro_amr_prolong, hydro_amr_prolong_action);

	void hydro_next_u(integer, real dt); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, hydro_next_u, hydro_next_u_action);

	std::pair<real,std::vector<real>> hydro_next_du(integer); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, hydro_next_du, hydro_next_du_action);

	void hydro_project(integer); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, hydro_project, hydro_project_action);

	//
	std::vector<real> hydro_get_bnd(integer, integer, exchange_type); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, hydro_get_bnd, hydro_get_bnd_action);
	//

	std::vector<real> hydro_get_amr_bnd(integer, integer, integer); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, hydro_get_amr_bnd, hydro_get_amr_bnd_action);
	//
	std::vector<real> hydro_restrict(integer); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, hydro_restrict, hydro_restrict_action);
	//

};

#endif

