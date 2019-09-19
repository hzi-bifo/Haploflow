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
#include <boost/range/combine.hpp>
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

// output operator for std::vector
template<class T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
    {
        os << "[ ";
        for (auto&& elem : v)
        {
            os << elem << " ";
        }
        os << "]";
        return os;
    }

struct Visits {
    std::vector<unsigned int> visits;
    friend std::ostream& operator<<(std::ostream& os, const Visits& v)
    {
        os << "[ ";
        for (auto&& elem : v.visits)
        {
            os << elem << " ";
        }
        os << "]";
        return os;
    };
};

struct EdgeProperties {
    std::string name;
    float capacity;
    float residual_capacity;
    Capacity cap_info;
    bool visited;
    bool first_vertex;
    std::vector<unsigned int> visits;
    Visits v; //DEBUG ONLY
    unsigned int last_visit;
    UEdge prev; // edge on path forwards
    UEdge next; // edge on path backwards
    float fatness; //fatness of the path going through edges forwards
    unsigned int distance; //distance from seed forwards
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
	UnitigGraph(deBruijnGraph&, std::string, float); // TODO delete dBg after UnitigGraph creation?
	UnitigGraph(); // debug
    ~UnitigGraph();
	void debug(); // debug information
    void assemble(std::string);
    void printGraph(std::ostream&, unsigned int cc);
    void dijkstra(UEdge seed, bool residual, bool local, unsigned int cc);
private:
	void connectUnbalanced(Vertex*, unsigned int*, std::string, deBruijnGraph&, float, float threshold);
	std::vector<std::pair<Vertex*,std::string> > addNeighbours(std::string& curr, const std::vector<char>&, const std::vector<char>&, deBruijnGraph&, unsigned int*, UVertex&, float threshold);
	std::pair<Vertex*,std::string> buildEdge(UVertex, Vertex*, std::string, std::string&, unsigned int*, float, float, deBruijnGraph&, float, float threshold);
	std::pair<Vertex*,std::string> buildEdgeReverse(UVertex, Vertex*, std::string, std::string&, unsigned int*, float, float, deBruijnGraph&, float, float threshold);
	UVertex addVertex(unsigned int*, std::string name, unsigned int ccc);
    std::vector<UEdge> get_sources(unsigned int cc);

    std::pair<std::string, float> calculate_contigs(std::vector<UEdge>&, std::vector<float>&, unsigned int cc);
    void reduce_flow(std::vector<UEdge>&, float, std::vector<float>&, std::set<unsigned int>&, unsigned int cc);
	std::vector<UEdge> find_fattest_path(UEdge seed, unsigned int cc);
    std::vector<UEdge> fixFlow(UEdge, unsigned int cc);
    
	//UEdge getSeed() const;
	std::vector<float> calculate_thresholds(deBruijnGraph&, std::string, float);
    std::vector<float> get_thresholds(std::vector<std::map<unsigned int, unsigned int>>& cov_distr, std::string);
    std::vector<float> rolling(std::vector<float>& in, unsigned int len);
    std::vector<float> cummin(std::vector<float>& in, unsigned int pos);

    std::vector<UEdge> blockPath(UEdge, unsigned int, unsigned int cc);
    std::vector<float> find_paths(unsigned int cc);
    std::pair<UEdge, bool> checkUnvisitedEdges(UEdge, unsigned int cc);
    std::pair<UEdge, bool> getUnvisitedEdge(const std::vector<UEdge>&, unsigned int, unsigned int cc);
    float remove_non_unique_paths(std::vector<std::vector<UEdge>>&, std::vector<UEdge>&, unsigned int, unsigned int, unsigned int cc);
    std::pair<UEdge, float> get_target(UEdge, bool, unsigned int cc);
    UEdge get_next_source(unsigned int cc);
    
    void cleanPath(std::vector<UEdge>&, unsigned int cc);
    void cleanGraph(unsigned int cc);
	void removeStableSets(unsigned int cc);
	void contractPaths(unsigned int cc);
    void removeEmpty(unsigned int cc);
    bool hasRelevance(unsigned int cc);
    void unvisit(unsigned int cc);

    float in_capacity(UVertex, unsigned int cc);
    float out_capacity(UVertex, unsigned int cc);

	bool test_hypothesis(float to_test_num, float to_test_denom, float h0, float threshold);

	unsigned int cc_; // used to mark the CC's. Since some of them might be deleted later on, does not represent the number of cc's
	std::vector<UGraph*> graphs_;
	std::vector<std::unordered_map<unsigned int, UVertex>> graph_map_;

    std::vector<float> thresholds_; // TODO
	
};

#endif
