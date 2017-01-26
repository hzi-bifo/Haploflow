#ifndef UNIG_H
#define UNIG_H

#include <boost/graph/graphviz.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/filtered_graph.hpp>
#include "deBruijnGraph.h"
#include <unordered_set>

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
						boost::property<boost::vertex_index1_t, unsigned int,
						boost::property<boost::vertex_name_t, std::string> >,
							boost::property<boost::edge_name_t, std::string,
							boost::property<boost::edge_capacity_t, float,
							boost::property<boost::edge_residual_capacity_t, unsigned int> > > > UGraph;

typedef typename boost::graph_traits<UGraph>::vertex_descriptor UVertex;
typedef typename boost::graph_traits<UGraph>::edge_descriptor UEdge;
typedef boost::graph_traits<UGraph>::vertex_iterator uvertex_iter;

// wrapper around the unitig graph defined above
class UnitigGraph {
public:
	// create a UnitigGraph from a dBg and its unbalanced vertices
	UnitigGraph(deBruijnGraph&);
private:
	unsigned int connectUnbalanced(Vertex*, unsigned int, std::string, deBruijnGraph&);
	unsigned int addNeighbours(std::string& curr, const std::vector<char>&, const std::vector<char>&, deBruijnGraph&, unsigned int, UVertex&);
	bool buildEdge(UVertex, Vertex*, std::string, std::string&, unsigned int, unsigned int, deBruijnGraph&);
	bool buildEdgeReverse(UVertex, Vertex*, std::string, std::string&, unsigned int, unsigned int, deBruijnGraph&);
	UVertex addVertex(unsigned int, std::string name);
	void cleanGraph();

	UGraph g_;
	std::unordered_map<unsigned int, UVertex> graph_;
};

#endif
