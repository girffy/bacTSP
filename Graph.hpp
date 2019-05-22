#ifndef GRAPH_HPP
#define GRAPH_HPP
#include "glpk.h"
#include "TSP_prob.hpp"
#include <vector>
#include <utility>

using namespace std;

/* a class for representing a sparse weighted, undirected graph via adjacency
 * lists, and providing some graph theory magic
 */
using algraph = vector<vector<pair<int, double>>>;

class Graph {
  public:
    // initialize a graph from the primal variables of tp.lp
    Graph(int num_nodes, TSP_prob tp);

    // initialize a graph from an adjacency list
    Graph(int num_nodes, vector<vector<pair<int, double>>> al);

    // compute the minimum-weight nontrivial cut using stoer-wagner
    // note: this destroys the graph!
    // TODO: this function assumes there is at least one nonzero edge; this is
    // valid? for TSP use cases, since every vertex has degree 2
    // if lim is supplied, if the algorithm finds a cut of weight at most lim,
    // it will return that cut instead of searching for a smaller one, thus
    // returning faster at the cost of cut optimality
    pair<vector<bool>, double> min_cut(double lim);

  private:
    // number of nodes and edges, respectively
    int n, m;

    // adjacency lists; G[i][j] == {n, w} indicates an edge (i,n) with weight w
    vector<vector<pair<int, double>>> G;



    // mutate this graph, merging vertices s and t. The vertex t is removed, the
    // edge {s,t} is removed if it exists, and for every other vertex v, the
    // weight of the edge {t,v} is added to the weight of {s,v}, creating the
    // new edge if needed.
    void merge(int s, int t);

    // recursive implementation of stoer-wagner; called by min_cut; nodes is the
    // set of nodes still under consideration
    pair<vector<bool>, double> min_cut_rec(double lim, vector<int> nodes);
};

#endif
