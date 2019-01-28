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

struct EdgeProperties;
struct VertexProperties;

struct Capacity {
    float first;
    float last;
    float min;
    float max;
    float avg;
    float length;
    float starting;
    float ending;
    friend std::ostream& operator<<(std::ostream& os, const Capacity& cap)
    {
        return os << "(" << cap.avg << ", " << cap.first << "/" << cap.last << ", (length: " << cap.length << "), (min: " << cap.min << ", max: " << cap.max << ")";
    }
};


typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
							VertexProperties,
							EdgeProperties> UGraph;

typedef typename boost::graph_traits<UGraph>::vertex_descriptor UVertex;
typedef typename boost::graph_traits<UGraph>::edge_descriptor UEdge;
typedef boost::graph_traits<UGraph>::vertex_iterator uvertex_iter;

struct Visits {
    std::vector<unsigned int> visits;
    friend std::ostream& operator<<(std::ostream& os, const Visits& visits)
    {
        os << "[ ";
        for (auto&& v : visits.visits)
        {
            os << v << " ";
        }
        os << "]";
        return os;
    }
};

struct EdgeProperties {
    std::string name;
    float capacity;
    float residual_capacity;
    Capacity cap_info;
    bool visited;
    //std::vector<unsigned int> visits;
    Visits visits;
    unsigned int last_visit;
    UEdge prev; // edge on path forwards
    UEdge next; // edge on path backwards
    float fatness; //fatness of the path going through edges forwards
    float fatness2; //fatness of the path going through edges backwards
    unsigned int distance; //distance from seed forwards
    unsigned int distance2; //distance from seed backwards
};

struct VertexProperties {
    unsigned int index;
    unsigned int scc;
    unsigned int cc;
    std::string name;
    unsigned int tarjan_index;
    bool onStack;
    int visiting_time;
};

typedef std::vector<UVertex> Connected_Component; // to distinguish from regular std::vector<UVertex>

// wrapper around the unitig graph defined above
class UnitigGraph {
public:
	// create a UnitigGraph from a dBg and its unbalanced vertices
	UnitigGraph(deBruijnGraph&, float); // TODO delete dBg after UnitigGraph creation?
	UnitigGraph(); // debug
	void debug(); // debug information
    void assemble(std::string);
    void printGraph(std::ostream&) const;
    void dijkstra(UEdge seed, bool);
private:
	void connectUnbalanced(Vertex*, unsigned int*, std::string, deBruijnGraph&, float);
	std::vector<std::pair<Vertex*,std::string> > addNeighbours(std::string& curr, const std::vector<char>&, const std::vector<char>&, deBruijnGraph&, unsigned int*, UVertex&);
	std::pair<Vertex*,std::string> buildEdge(UVertex, Vertex*, std::string, std::string&, unsigned int*, float, float, deBruijnGraph&, float);
	std::pair<Vertex*,std::string> buildEdgeReverse(UVertex, Vertex*, std::string, std::string&, unsigned int*, float, float, deBruijnGraph&, float);
	UVertex addVertex(unsigned int*, std::string name);

    std::pair<std::string, std::pair<float, float> > calculate_contigs(std::vector<UEdge>&);
	float calculate_thresholds(const deBruijnGraph&, float);
    float calculate_gain(UVertex& v);
    std::pair<float, std::vector<float> > calculate_flow(std::vector<UEdge>&);
	std::vector<UEdge> find_fattest_path(UEdge);
    
	UEdge getSeed() const;

    void blockPath(UEdge, unsigned int);
    void fixFlow();
    std::string contig_from_blocked_path(UEdge, unsigned int);
    std::pair<UEdge, bool> checkUnvisitedEdges(UEdge);
    std::pair<UEdge, bool> getUnvisitedEdge(const std::vector<UEdge>&, unsigned int);
    void markCycles();
    void cleanGraph();
	void removeStableSets();
	void contractPaths();
    void removeEmpty();
    void unvisit();

    float in_capacity(UVertex);
    float out_capacity(UVertex);

	bool test_hypothesis(float to_test_num, float to_test_denom, float h0);

	unsigned int cc_; // used to mark the CC's. Since some of them might be deleted later on, does not represent the number of cc's
	UGraph g_;
	std::unordered_map<unsigned int, UVertex> graph_;
    float threshold_; // stores a treshold for edges to consider
	
};

#endif
