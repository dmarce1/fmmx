/*
 * silo_output.hpp
 *
 *  Created on: May 31, 2014
 *      Author: dmarce1
 */

#ifndef SILO_OUTPsUT_HPP_
#define SILO_OUTPsUT_HPP_

#ifndef NO_OUTPUT
#include <silo.h>
#endif

#include "defs.hpp"

class node_client;

#include <hpx/lcos/local/mutex.hpp>
#include <vector>
#include <hpx/runtime/actions/component_action.hpp>
#include <hpx/runtime/components/server/managed_component_base.hpp>

class silo_output: public hpx::components::managed_component_base<silo_output> {
public:
	static constexpr double precision = 1.0e-10;
	static constexpr int Nchild = 1 << NDIM;
	struct zone {
		std::array<real, NF> fields;
		std::array<real, NDIM> position;
		std::array<real, NDIM> span;
		zone() {
		}
		zone(const zone& z) {
			*this = z;
		}
		zone(zone&& z) {
			*this = z;
		}
		zone& operator=(const zone& z) {
			fields = z.fields;
			position = z.position;
			span = z.span;
			return *this;
		}
		zone& operator=(zone&& z) {
			fields = std::move(z.fields);
			position = std::move(z.position);
			span = std::move(z.span);
			return *this;
		}
		template<typename Archive>
		void serialize(Archive& ar, const int v) {
			ar & fields;
			ar & position;
			ar & span;
		}
	};
	struct silo_zone {
		std::array<real,NF> fields;
		std::vector<int> vertices;
		silo_zone() :
				vertices(Nchild) {
		}
		silo_zone(const silo_zone& s) :
				vertices(Nchild) {
			fields = s.fields;
			vertices = s.vertices;
		}
		silo_zone(silo_zone&& s) :
				vertices(Nchild) {
			fields = std::move(s.fields);
			vertices = std::move(s.vertices);
		}
	};
	struct vertex: public std::vector<double> {
		int index;
		vertex() :
				std::vector<double>(NDIM) {
		}
		vertex(const vertex& v) :
				std::vector<double>(v) {
			index = v.index;
		}
		vertex(vertex&& v) :
				std::vector<double>(NDIM) {
			std::vector<double>::operator=(std::move(*(static_cast<std::vector<double>*>(&v))));
			index = v.index;
		}
	};
	struct vertex_less_functor: std::binary_function<vertex, vertex, bool> {
		bool operator()(const vertex& x, const vertex& y) const {
			for (int i = 0; i < NDIM; i++) {
				if (x[i] - y[i] > precision) {
					return false;
				} else if (x[i] - y[i] < -precision) {
					return true;
				}
			}
			return false;
		}
	};
	using vertex_dir_type = std::set<vertex,vertex_less_functor>;
	using silo_zone_dir_type = std::list<silo_zone>;
private:
	std::vector<std::string> names;
	int current_index;
	vertex_dir_type nodedir;
	silo_zone_dir_type zonedir;
	mutable hpx::lcos::local::spinlock mutex0;
public:
	silo_output() = default;
	virtual ~silo_output() = default;
	void do_output(std::list<std::size_t> leaves);HPX_DEFINE_COMPONENT_ACTION_TPL( silo_output, do_output, do_output_action );
}
;

#endif /* SILO_OUTPUT_HPP_ */
