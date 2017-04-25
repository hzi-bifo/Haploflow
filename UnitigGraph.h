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
#include <algorithm> // heap
#include <queue>

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
						boost::property<boost::vertex_index1_t, unsigned int, // index
						boost::property<boost::vertex_index2_t, unsigned int, // cc
						boost::property<boost::vertex_name_t, std::string, // kmer
						boost::property<boost::vertex_discover_time_t, bool> > > >, // visited flag
							boost::property<boost::edge_name_t, std::string, // edge sequence
							boost::property<boost::edge_index_t, bool, //edge visited flag
							boost::property<boost::edge_capacity_t, float, // capacity
							boost::property<boost::edge_residual_capacity_t, float> > > > > UGraph;

typedef typename boost::graph_traits<UGraph>::vertex_descriptor UVertex;
typedef typename boost::graph_traits<UGraph>::edge_descriptor UEdge;
typedef boost::graph_traits<UGraph>::vertex_iterator uvertex_iter;

typedef std::vector<UVertex> Connected_Component; // to distinguish from regular std::vector<UVertex>

// wrapper around the unitig graph defined above
class UnitigGraph {
public:
	// create a UnitigGraph from a dBg and its unbalanced vertices
	UnitigGraph(deBruijnGraph&); // TODO delete dBg after UnititgGraph creation?
	void calculateFlow();
private:
	void connectUnbalanced(Vertex*, unsigned int*, std::string, deBruijnGraph&);
	std::vector<std::pair<Vertex*,std::string> > addNeighbours(std::string& curr, const std::vector<char>&, const std::vector<char>&, deBruijnGraph&, unsigned int*, UVertex&);
	std::pair<Vertex*,std::string> buildEdge(UVertex, Vertex*, std::string, std::string&, unsigned int*, float, float, deBruijnGraph&);
	std::pair<Vertex*,std::string> buildEdgeReverse(UVertex, Vertex*, std::string, std::string&, unsigned int*, float, float, deBruijnGraph&);
	UVertex addVertex(unsigned int*, std::string name);
	void cleanGraph();
	void contractPaths();
	
	void find_fattest_path(UVertex target, std::string& sequence, std::vector<float>& coverage_fraction, std::vector<UEdge>& visited_edges);
	void add_sorted_edges(std::vector<UEdge>& q, const UVertex& source, bool addAll);
	bool test_hypothesis(float to_test, float h0);
	std::vector<std::vector<UVertex> > getSources() const;

	unsigned int cc_; // used to mark the CC's. Since some of them might be deleted later on, does not represent the number of cc's
	UGraph g_;
	std::unordered_map<unsigned int, UVertex> graph_;
	
};

#endif
