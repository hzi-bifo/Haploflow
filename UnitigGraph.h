#ifndef UNIG_H
#define UNIG_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/filtered_graph.hpp>
#include "deBruijnGraph.h"
#include <unordered_set>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
						boost::property<boost::vertex_name_t, std::string>,
							boost::property<boost::edge_name_t, std::string,
							boost::property<boost::edge_capacity_t, unsigned int,
							boost::property<boost::edge_residual_capacity_t, unsigned int> > > > UGraph;

typedef typename boost::graph_traits<UGraph>::vertex_descriptor UVertex;
typedef typename boost::graph_traits<UGraph>::edge_descriptor UEdge;

// wrapper around the unitig graph defined above
class UnitigGraph {
public:
	// create a UnitigGraph from a dBg and its unbalanced vertices
	UnitigGraph(const deBruijnGraph&);
private:
	UGraph g_;
};

#endif
