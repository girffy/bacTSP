#include "Graph.hpp"
#include <algorithm>
#include <vector>
#include <iostream>

#ifndef DEBUG
#define printf(...) (void)0
#endif

const double INF = numeric_limits<double>::infinity();
double EPS = 1e-10;

Graph::Graph(int num_nodes, TSP_prob tp){
  n = num_nodes;
  m = n*(n-1)/2;

  G.resize(n); // TODO: does this initialize the elements?
  for(int i = 1; i <= m; i++){
    auto [src,dest] = tp.var_to_edge(i);
    double x = glp_get_col_prim(tp.lp, i);
    if(x > EPS){
      G[src].push_back({dest, x});
      G[dest].push_back({src, x});
    }
  }
}

Graph::Graph(int num_nodes, vector<vector<pair<int, double>>> al){
  n = num_nodes;
  m = n*(n-1)/2;
  G = al;
}

// TODO: maybe there's a faster way to do this by linking the two "copies" of
// each edge together somehow? but .erase() is still expensive
void Graph::merge(int s, int t){
  //printf("Merging %d,%d\n", s, t);
  assert(s<n && t<n && s != t);

  // TODO: static
  // TODO: these are 0-initialized, right?
  vector<double> tweights(n);
  for(auto [nd, wt] : G[t]){
    tweights[nd] += wt;
  }


  for(int srcnode = 0; srcnode < (int)G.size(); srcnode++){
    if(srcnode == t) continue;
    bool found_t = false;
    bool found_s = false;
    for(int i = G[srcnode].size()-1; i >= 0; i--){
      int destnode = G[srcnode][i].first;
      if(destnode == t){
        found_t = true;
        G[srcnode].erase(G[srcnode].begin()+i, G[srcnode].begin()+i+1);
        if(i != s) continue;
      }else if(destnode == s){
        found_s = true;
        G[srcnode][i].second += tweights[srcnode];
      }else if(srcnode == s){
        G[srcnode][i].second += tweights[destnode];
      }
    }

    // if this node had an edge to t but not s, we need to add a new edge to s
    // with the appropriate weight
    if(srcnode != s && found_t && !found_s){
      G[srcnode].push_back({s, tweights[srcnode]});
      G[s].push_back({srcnode, tweights[srcnode]});
    }
  }

  G[t].clear();
}

void print_cut(vector<bool> cut){
  for(auto it = cut.begin(); it != cut.end(); it++){
    if(*it) printf("%d ", (int)(it-cut.begin()));
  }
}

void print_vlist(vector<int> A){
  for(int i : A) printf("%d ", i);
}

void print_flist(vector<double> xs){
  for(double x : xs){
    printf("%.2f ", x);
  }
}

void print_edges(vector<vector<pair<int, double>>> G){
  for(int i = 0; i < (int)G.size(); i++){
    printf("%d: ", i);
    for(auto [nd, wt] : G[i]) printf("%d(%.2f) ", nd, wt);
    printf("\n");
  }
}

pair<vector<bool>, double> Graph::min_cut(double lim){
  vector<int> nodes(n);
  for(int i = 0; i < n; i++) nodes[i] = i;
  //print_edges(G);
  return min_cut_rec(lim, nodes);
}

pair<vector<bool>, double> Graph::min_cut_rec(double lim, vector<int> nodes){
  //printf("ns=%ld\n", nodes.size());
  assert(nodes.size() >= 2);
  // TODO: choose randomly? or via some heuristic?
  vector<int> A = {nodes[0]};

  // TODO: static
  vector<double> scores(n); // TODO: 0-init?
  while(A.size() < nodes.size()){
    int last_addition = A[A.size()-1];
    scores[last_addition] = -INF;
    for(auto [nd, wt] : G[last_addition]){
      scores[nd] += wt;
    }
    //printf("scores: "); print_flist(scores); printf("\n");

    // TODO: make this more efficient with a heap of some kind
    //int next = std::min_element(scores.begin(), scores.end()) - scores.begin();
    int next = -1;
    double highest_score = -INF;
    for(int nd : nodes){
      if(scores[nd] > highest_score){
        next = nd;
        highest_score = scores[nd];
      }
    }
    assert(next != -1);
    A.push_back(next);
  }

  int s = A[A.size() - 2];
  int t = A[A.size() - 1];
  double cut_size = scores[t];
  //printf("A(%.2f) is: ", cut_size); print_vlist(A); printf("\n");
  if(nodes.size() > 2 && cut_size >= lim){
    merge(s,t);
    std::remove(nodes.begin(), nodes.end(), t);
    nodes.resize(nodes.size()-1);
    //printf("Finished merging\n");
    //print_edges(G);
    auto [rec_cut, rec_wt] = min_cut_rec(lim, nodes);
    if(rec_wt < cut_size){
      rec_cut[t] = rec_cut[s];
      //printf("ns=%ld; returning rec cut (", nodes.size()); print_cut(rec_cut); printf("); wt=%f\n", rec_wt);
      return {rec_cut, rec_wt};
    }
  }

  vector<bool> cut(n);
  cut[t] = true;
  //printf("ns=%ld; returning this cut (", nodes.size()); print_cut(cut); printf("); wt=%f\n", cut_size);
  return {cut, cut_size};
}

