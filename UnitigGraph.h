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

struct VertexProperties {
    unsigned int index;
    unsigned int scc;
    unsigned int cc;
    std::string name;
    unsigned int tarjan_index;
    bool onStack;
    bool visited;
};

struct Capacity {
    float first;
    float last;
    float min;
    float max;
    float avg;
    float length;
    friend std::ostream& operator<<(std::ostream& os, const Capacity& cap)
    {
        return os << "(" << cap.avg << ", " << cap.first << "/" << cap.last << ", (" << cap.length << "), (" << cap.min << ", " << cap.max << "))";
    }
};

struct EdgeProperties{
    std::string name;
    float capacity;
    Capacity cap_info;
    float residual_capacity;
    bool visited;
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
							VertexProperties,
							EdgeProperties> UGraph;

typedef typename boost::graph_traits<UGraph>::vertex_descriptor UVertex;
typedef typename boost::graph_traits<UGraph>::edge_descriptor UEdge;
typedef boost::graph_traits<UGraph>::vertex_iterator uvertex_iter;

typedef std::vector<UVertex> Connected_Component; // to distinguish from regular std::vector<UVertex>

// wrapper around the unitig graph defined above
class UnitigGraph {
public:
	// create a UnitigGraph from a dBg and its unbalanced vertices
	UnitigGraph(deBruijnGraph&, float); // TODO delete dBg after UnitigGraph creation?
	UnitigGraph(); // debug
	void debug(); // debug information
    void calculateFlow();
    void printGraph(std::ostream&) const;
private:
	void connectUnbalanced(Vertex*, unsigned int*, std::string, deBruijnGraph&, float);
	std::vector<std::pair<Vertex*,std::string> > addNeighbours(std::string& curr, const std::vector<char>&, const std::vector<char>&, deBruijnGraph&, unsigned int*, UVertex&);
	std::pair<Vertex*,std::string> buildEdge(UVertex, Vertex*, std::string, std::string&, unsigned int*, float, float, deBruijnGraph&);
	std::pair<Vertex*,std::string> buildEdgeReverse(UVertex, Vertex*, std::string, std::string&, unsigned int*, float, float, deBruijnGraph&);
	UVertex addVertex(unsigned int*, std::string name);

	float calculateThresholds(const deBruijnGraph&, float);
    void markCycles();
    void cleanGraph();
	void removeStableSets();
	void contractPaths();
    void removeEmpty();

    float in_capacity(UVertex);
    float out_capacity(UVertex);

	std::vector<UEdge> find_fattest_path(UEdge);
    std::pair<std::string, std::pair<float, float> > calculate_contigs(std::vector<UEdge>&);
	bool test_hypothesis(float to_test_num, float to_test_denom, float h0, bool = false);
    float calculate_gain(UVertex& v);
	UEdge getSeed() const;

	unsigned int cc_; // used to mark the CC's. Since some of them might be deleted later on, does not represent the number of cc's
	UGraph g_;
	std::unordered_map<unsigned int, UVertex> graph_;
    float threshold_; // stores a treshold for edges to consider
	
};

#endif
