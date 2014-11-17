/*
 * node_client.hpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dmarce1
 */

#ifndef NODE_CLIENT_HPP_
#define NODE_CLIENT_HPP_

#include <hpx/include/components.hpp>
#include <boost/serialization/vector.hpp>
#include "defs.hpp"

class node_client {
private:
	hpx::id_type id;
public:
	operator hpx::id_type() const;
	node_client();
	node_client(const hpx::id_type&);
	node_client& operator=(const hpx::id_type&);
	void set_boundary(hpx::future<std::vector<real>>& f, integer d);
	void set_multipoles(hpx::future<std::vector<real>>& f, integer ci);
	void set_expansions(hpx::future<std::vector<real>>&);
	hpx::future<std::vector<real>> get_data() const;
	hpx::future<integer> get_node_count() const;
	hpx::future<std::list<std::size_t>> get_leaf_list() const;
	template<class Arc>
	void serialize(Arc&, const unsigned);
};

template<class Arc>
void node_client::serialize(Arc& a, const unsigned v) {
	boost::serialization::serialize(a, id, v);
}

#endif /* NODE_CLIENT_HPP_ */
