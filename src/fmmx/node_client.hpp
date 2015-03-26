/*
 * node_client.hpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#ifndef NODE_CLIENT_HPP_
#define NODE_CLIENT_HPP_

#include "defs.hpp"
#include <boost/serialization/vector.hpp>

class node_client {
private:
	hpx::id_type id;
public:
	hpx::future<void> get_tree(std::vector<node_client> my_neighbors = std::vector<node_client>(NNEIGHBOR, node_client(hpx::invalid_id)));
	hpx::future<std::vector<real>> hydro_restrict(integer rk);
	hpx::future<bool> refine();
	hpx::future<void> refine_proper();
	hpx::future<void> set_me(hpx::id_type id);
	hpx::future<void> derefine();
	hpx::future<void> destroy();
	operator hpx::id_type() const;
	node_client();
	hpx::future<void> hydro_exchange(integer); //
	hpx::future<void> hydro_amr_prolong(integer); //
	hpx::future<std::vector<real>> hydro_get_amr_bnd(integer, integer, integer); //
	hpx::future<std::vector<real>> hydro_get_bnd(integer, integer); //
	hpx::future<std::vector<node_client>> get_children() const;
	node_client(const hpx::id_type&);
	node_client& operator=(const hpx::id_type&);
	bool operator==(const hpx::id_type&) const;
	bool operator!=(const hpx::id_type&) const;
	void set_fmm_boundary(hpx::future<std::vector<real>>& f, integer d);
	void set_multipoles(hpx::future<std::vector<real>>& f, integer ci);
	void set_expansions(hpx::future<std::vector<real>>&);
	hpx::future<std::vector<double>> get_data() const;
	hpx::future<integer> get_node_count() const;
	hpx::future<std::list<std::size_t>> get_leaf_list() const;
	hpx::future<real> execute(integer rk);
	template<class Arc>
	void serialize(Arc&, const unsigned);

	hpx::future<void> hydro_next_u(integer, real dt); //
	hpx::future<std::pair<real, std::vector<real>>> hydro_next_du(integer); //
	hpx::future<void> hydro_project(integer); //

};

template<class Arc>
void node_client::serialize(Arc& a, const unsigned v) {
	boost::serialization::serialize(a, id, v);
}

#endif /* NODE_CLIENT_HPP_ */
