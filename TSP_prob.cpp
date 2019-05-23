#include "TSP_prob.hpp"
#include "Graph.hpp"

#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>
#include <regex>

#ifndef DEBUG
#define printf(...) (void)0
#endif
#define logf printf("%*c", depth*2, ' '); printf

TSP_prob::TSP_prob(string filename) {
  ifstream infile(filename);
  assert(!infile.fail());
  string line;
  regex prop_regex("^([^\\s:]+)\\s*:\\s+(.*)$");
  map<string, string> props;

  // read metadata
  while(getline(infile, line)){
    smatch match;
    bool matched = regex_match(line, match, prop_regex);
    if(!matched) break;

    props.insert({match[1], match[2]});
    infile.clear();
  }

  num_nodes = stoi(props["DIMENSION"]);
  edges = new edge[num_nodes*(num_nodes-1)/2];

  // read node/edge data into this->d
  assert(line == "NODE_COORD_SECTION");
  if(props["EDGE_WEIGHT_TYPE"] == "EUC_2D"){
    vector<pair<double, double>> nodes(num_nodes);
    while(getline(infile, line)){
      if(line == "EOF") break;
      istringstream iss(line);
      int idx;
      double x, y;
      iss >> idx >> x >> y;

      nodes[idx-1] = {x,y};
    }


    // populate distance matrix
    this->d.resize(num_nodes);
    for(int i = 0; i < num_nodes; i++){
      this->d[i] = vector<double>(num_nodes);
      for(int j = 0; j < num_nodes; j++){
        double dx = nodes[i].first - nodes[j].first;
        double dy = nodes[i].second - nodes[j].second;

        // per TSPLIB standard, distances are rounded to integers
        d[i][j] = round(sqrt(dx*dx + dy*dy));
      }
    }
  }
  init_LP();
}

TSP_prob::TSP_prob(dmtx d) : d(d) {
  int n = d.size();
  this->num_nodes = n;
  this->edges = new edge[n*(n-1)/2];
  this->init_LP();
}

TSP_prob::~TSP_prob(){
  delete[] this->edges;
  //printf("Deleting %p\n", (void*)this->lp);
  if(this->lp) {
      glp_delete_prob(this->lp);
  }
}

TSP_prob::TSP_prob(const TSP_prob &other){
  this->lp = glp_create_prob();
  //printf("Creating1 %p\n", (void*)this->lp);
  glp_copy_prob(this->lp, other.lp, GLP_ON);

  int n = other.num_nodes;
  int ne = n*(n-1)/2;

  this->num_nodes = n;
  this->d = other.d;
  this->lb = other.lb;
  this->ub = other.ub;
  this->depth = other.depth;
  this->edges = new edge[ne];
  std::copy(other.edges, other.edges+ne, this->edges);
}

TSP_prob& TSP_prob::operator=(TSP_prob&& other){
  this->num_nodes = other.num_nodes;
  this->d = other.d;
  this->lp = other.lp;
  this->lb = other.lb;
  this->ub = other.ub;
  this->depth = other.depth;
  this->edges = other.edges;

  other.lp = NULL;
  other.edges = NULL;

  return *this;
}

void TSP_prob::init_LP(){
  int n = this->d.size();
  auto ia = new int[1 + n*n];
  auto ja = new int[1 + n*n];
  auto ar = new double[1 + n*n];
  lp = glp_create_prob();
  //printf("Creating2 %p\n", (void*)lp);
  glp_set_prob_name(lp, "TSP");
  glp_set_obj_dir(lp, GLP_MIN);
  glp_add_rows(lp, n);
  glp_add_cols(lp, n*(n-1)/2);

  // add rows/constraints
  constexpr auto BUFSZ = 50;
  char buffer[BUFSZ];
  for(int i=0; i<n; i++){
    if (auto sz = snprintf(buffer, BUFSZ, "degcon%d", i);
        sz < 0 || sz >= BUFSZ) { exit (-2); }
    glp_set_row_name(lp, i+1, buffer);
    glp_set_row_bnds(lp, i+1, GLP_FX, 2.0, 2.0);
  }

  // add variables
  int vidx = 1;
  for(int i=0; i<n; i++){
    for(int j=i+1; j<n; j++){
      if (auto sz = snprintf(buffer, BUFSZ, "x%d,%d", i, j);
          sz < 0 || sz >= BUFSZ) { exit (-2); }
      glp_set_col_name(lp, vidx, buffer);
      glp_set_col_bnds(lp, vidx, GLP_DB, 0.0, 1.0);
      glp_set_obj_coef(lp, vidx, d[i][j]);
      this->edges[vidx-1] = {i,j};
      vidx++;
    }
  }

  // populate constraint matrix
  int aidx = 1;
  vidx = 1;
  for(int a=0; a<n; a++){
    for(int b=a+1; b<n; b++){
      ia[aidx] = a+1;
      ja[aidx] = vidx;
      ar[aidx] = 1.0;
      aidx++;

      ia[aidx] = b+1;
      ja[aidx] = vidx;
      ar[aidx] = 1.0;
      aidx++;
      vidx++;
    }
  }

  glp_load_matrix(lp, aidx-1, ia, ja, ar);
  this->lp = lp;

  delete[] ia;
  delete[] ja;
  delete[] ar;
}

void TSP_prob::add_subtour_constraint(vector<bool> cut){
  logf("Adding subtour constraint:");
  if(num_nodes < 50){
    for(int i = 0; i < num_nodes; i++){
      if(cut[i]) printf(" %d", i);
    }
  }
  printf("\n");

  int n = this->num_nodes;
  // TODO: make ind and val scratch variables in TSP_prob
  int *ind = new int[1 + n*n];
  double *val = new double[1 + n*n];
  int eidx = 1;
  glp_add_rows(this->lp, 1);
  int r = glp_get_num_rows(this->lp);
  glp_set_row_name(this->lp, r, "subtour");
  glp_set_row_bnds(this->lp, r, GLP_LO, 2.0, 0.0);

  for(int src = 0; src < n; src++){
    for(int dest = src+1; dest < n; dest++){
      if(cut[src] != cut[dest]){
        int idx = edge_to_var({src,dest});
        ind[eidx] = idx;
        val[eidx] = 1.0;
        eidx++;
      }
    }
  }

  glp_set_mat_row(this->lp, r, eidx-1, ind, val);

  delete[] ind;
  delete[] val;
}

void TSP_prob::fix_var(int var, double value){
  int ind[1+1];
  double val[1+1];
  char buffer[50];
  auto [src, dest] = var_to_edge(var);

  glp_add_rows(this->lp, 1);
  int r = glp_get_num_rows(this->lp);
  sprintf(buffer, "x%d,%d=%f", src, dest, value);
  glp_set_row_name(this->lp, r, buffer);
  glp_set_row_bnds(this->lp, r, GLP_FX, value, value);

  ind[1] = var;
  val[1] = 1.0;
  glp_set_mat_row(this->lp, r, 1, ind, val);
}

int TSP_prob::edge_to_var(edge e){
  auto [src, dest] = e;
  return this->edge_to_var(src, dest);
}

// This function in broken in ways I dont understand how to fix (columns shouldnt be negative idx)
// do it the naive way until b2 fixes the tricky way
int TSP_prob::edge_to_var(int src, int dest){
  assert(src < dest);
  // idx is diff of partial sums + offs
  auto psum = [](auto n) { return n*(n-1)/2; };
  int idx = psum(this->num_nodes) - psum(this->num_nodes-src) + (dest-src);
  // printf("EDGE (%d, %d) TO VAR %d\n", src, dest, idx);
  return idx;
}

edge TSP_prob::var_to_edge(int i){
  return this->edges[i-1];
}


const double EPS = 1e-10;
const double INF = numeric_limits<double>::infinity();

void TSP_prob::solve(){
  this->lb = 0.0;
  this->ub = INF;
  this->depth = 0;
  bool success = this->recsolve();
  assert(success);
}

bool TSP_prob::recsolve(){
  logf("recsolve: depth=%d, range=[%.3f,%.3f]\n", depth, lb, ub);
  int n = this->num_nodes;
  bool feas = this->cp_solve();
  double obj_val = glp_get_obj_val(lp);
  if(!feas || glp_get_obj_val(lp) > ub || lb > ub-EPS){
    if(!feas){
      logf("fail: cp_solve could not find feas\n");
    }else if(glp_get_obj_val(lp) > ub){
      logf("fail: tour length too high (%f>%f)\n", glp_get_obj_val(lp), ub);
    }else if(lb > ub-EPS){
      logf("fail: bound range is empty [%f, %f]\n", lb, ub);
    }
    return false;
  }
  this->lb = max(this->lb, obj_val);

  // find the most fractional primal variable, or determine that the solution is
  // integral
  int branch_var = -1;
  double record = 10.0;
  for(int i = 0; i < n*(n-1)/2; i++){
    double x = glp_get_col_prim(lp, i+1);
    if(EPS < x && x < 1-EPS && (branch_var == -1 || abs(0.5-x) < record)){
      branch_var = i+1;
      record = abs(0.5-x);
    }
  }
  if(branch_var == -1){
    ub = glp_get_obj_val(lp);
    logf("found integer soln; ov=%.2f\n", glp_get_obj_val(lp));
    return true;
  }
  auto [src, dest] = var_to_edge(branch_var);

  // now, branch on branch_var
  // TODO: could be slightly more efficient by using the current TSP_prob for
  // one of the branches, instead of making two new copies

  // x=1 branch
  TSP_prob tp1(*this);
  tp1.depth++;
  tp1.fix_var(branch_var, 1.0);
  logf("branching on x%d,%d = 1\n", src, dest);
  bool feas1 = tp1.recsolve();
  double objval1 = feas1 ? glp_get_obj_val(tp1.lp) : INF;
  if(feas1){
    logf("ub: %.3f -> ", ub);
    this->ub = min(this->ub, tp1.ub);
    printf("%.3f\n", ub);
  }

  // x=0 branch
  auto& tp0 = *this;
  tp0.depth++;
  tp0.fix_var(branch_var, 0.0);
  logf("branching on x%d,%d = 0\n", src, dest);
  bool feas0 = tp0.recsolve();
  double objval0 = feas0 ? glp_get_obj_val(tp0.lp) : INF;

  logf("objvals are %f, %f\n", objval0, objval1);

  if(objval0 < objval1){
    logf("choosing branch x%d,%d = 0\n", src, dest);
    //*this = std::move(tp0);
  }else if(feas1){
    logf("choosing branch x%d,%d = 1\n", src, dest);
    *this = std::move(tp1);
  }else{
    logf("both branches failed\n");
    return false;
  }

  return true;
}

bool TSP_prob::cp_solve(){
  // solve and add subtour inequalities until none are violated
  logf("cp_solve enter\n");
  while(true){
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    assert(glp_simplex(this->lp, &parm) == 0);
    int status = glp_get_status(this->lp);
    if(status != GLP_OPT){
      assert(status == GLP_NOFEAS);
      logf("NOFEAS\n");
      return false;
    }

#ifdef BOOSTCUT
    auto [cut, weight] = this->min_cut();
#else
    Graph G(num_nodes, *this);
    auto [cut, weight] = G.min_cut(2.0 - EPS);
    //auto [cut, weight] = G.min_cut(0);
#endif
    //printf("weight is %.2f\n", weight);
    if(weight >= 2-EPS) break;
    this->add_subtour_constraint(cut);
    //print_wmat();
  }

  // TODO: look for more cutting planes, e.g. comb inequalities

  logf("cp_solve: found soln with ov=%f\n", glp_get_obj_val(lp));
  return true;
}


// debug function; print out the current weights in matrix form
void TSP_prob::print_wmat(){
  printf("    ");
  for(int k = 0; k < num_nodes; k++) printf("%-4d", k);
  printf("\n");
  for(int i = 0; i < num_nodes; i++){
    printf("%-4d", i);
    for(int j = 0; j < num_nodes; j++){
      if(i == j){
        printf("___ ");
      }else{
        edge e = {min(i,j), max(i,j)};
        double x = glp_get_col_prim(lp, edge_to_var(e));
        if(x != 0){
          printf("%.1f ", x);
        }else{
          printf("    ");
        }
      }
    }
    printf("\n");
  }
}

#ifdef BOOSTCUT
pair<vector<bool>, double> TSP_prob::min_cut(){
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property<boost::edge_weight_t, double>> undirected_graph;
  typedef boost::graph_traits<undirected_graph>::vertex_descriptor vertex_descriptor;
  typedef boost::property_map<undirected_graph, boost::edge_weight_t>::type weight_map_type;
  typedef boost::property_traits<weight_map_type>::value_type weight_type;

  int ne = this->num_nodes*(this->num_nodes-1)/2;
  double *wbuf = new double[ne]; // TODO: put this back to being a class variable

  for(int i = 0; i < ne; i++){
    wbuf[i] = glp_get_col_prim(this->lp, i+1);
  }


  undirected_graph g(edges, edges + ne, wbuf, this->num_nodes, ne);
  BOOST_AUTO(parities, boost::make_one_bit_color_map(num_vertices(g), get(boost::vertex_index, g)));
  double w = boost::stoer_wagner_min_cut(g, get(boost::edge_weight, g), boost::parity_map(parities));

  delete[] wbuf;

  vector<bool> v(num_nodes);
  for(int i = 0; i < num_nodes; i++){
    v[i] = get(parities, i);
  }

  return {v, w};
}
#endif

// could be made faster, but not exactly perf critical anyway
tour TSP_prob::get_tour(){
  assert(num_nodes > 2);
  tour t = {0};

  for(int i = 0; i < num_nodes-1; i++){
    for(int j = 0; j < num_nodes; j++){
      int cur = t[t.size()-1];
      if(cur == j) continue;
      edge e = {min(cur,j), max(cur,j)};
      int idx = edge_to_var(e);
      double x = glp_get_col_prim(lp, idx);
      assert(x < EPS || x > 1-EPS);
      if(x > 1-EPS && (t.size() == 1 || t[t.size()-2] != j)){
        t.push_back(j);
        break;
      }
    }
  }

  assert(t.size() == (unsigned int)num_nodes);

  return t;
}

void TSP_prob::dump_LP_soln(){
  printf("LP dump:\n");
  for(int eidx = 1; eidx <= glp_get_num_cols(lp); eidx++){
    edge e = var_to_edge(eidx);
    printf("\tx%d,%d = %f\n", e.first, e.second, glp_get_col_prim(lp, eidx));
  }
}
