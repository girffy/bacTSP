#include "TSP_prob.hpp"
#include <assert.h>
#include <stdio.h>
#include <algorithm>

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
  this->edges = new edge[ne];
  std::copy(other.edges, other.edges+ne, this->edges);
}

TSP_prob& TSP_prob::operator=(TSP_prob&& other){
  this->num_nodes = other.num_nodes;
  this->d = other.d;
  this->lp = other.lp;
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

void TSP_prob::add_subtour_constraint(parity_map pmap){
  printf("Adding subtour constraint:");
  for(int i = 0; i < num_nodes; i++){
    if(get(pmap, i)) printf(" %d", i);
  }
  printf("\n");

  int n = this->num_nodes;
  // TODO: make ind and val scratch variables in TSP_prob
  int *ind = new int[1 + n*n];
  double *val = new double[1 + n*n];
  int eidx = 1;
  glp_add_rows(this->lp, 1);
  int r = glp_get_num_rows(this->lp);
  printf("r is %d\n", r);
  glp_set_row_name(this->lp, r, "subtour");
  glp_set_row_bnds(this->lp, r, GLP_LO, 2.0, 0.0);

  double debug_weight = 0;
  for(int src = 0; src < n; src++){
    for(int dest = src+1; dest < n; dest++){
      if(get(pmap, src) != get(pmap, dest)){
        int idx = edge_to_var({src,dest});
        debug_weight += glp_get_col_prim(lp, idx);
        ind[eidx] = idx;
        val[eidx] = 1.0;
        eidx++;
      }
    }
  }

  glp_set_mat_row(this->lp, r, eidx-1, ind, val);

  /* big debug block */
  printf("debug_weight=%f\n", debug_weight);
  printf("eidx=%d\n", eidx);
  int rl = glp_get_mat_row(lp, r, ind, val);
  printf("rl is %d\n", rl);
  double *row = new double[1+n*n];
  for(int i = 0; i < 1+n*n; i++) row[i] = 0;
  for(int i = 1; i <= rl; i++){
    row[ind[i]] = val[i];
  }
  /*
  printf("new row is:\n");
  for(int i = 1; i <= glp_get_num_cols(lp); i++){
    auto [src, dest] = var_to_edge(i);
    //printf("%.2f ", row[i]);
    //printf("%d ", (int)row[i]);
    printf("%d,%d: %.2f\n", src, dest, row[i]);
  }
  printf("\n\n");
  */
  /*
  printf("current soln:\n");
  dump_LP_soln();
  printf("@@@@@@@@@@@@@@@@@@@@\n");
  printf("2nd last row primal: %f\n", glp_get_row_dual(lp, glp_get_num_rows(lp)-2));
  printf("2nd last row dual: %f\n", glp_get_row_dual(lp, glp_get_num_rows(lp)-2));
  printf("@@@@@@@@@@@@@@@@@@@@\n");
  */

  delete[] row;

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
  glp_prob *lp = this->lp;
  int n = this->num_nodes;
  bool feas = this->cp_solve();
  double obj_val = glp_get_obj_val(lp);
  if(!feas || glp_get_obj_val(lp) > ub || this->lb > ub-EPS){
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
  if(branch_var == -1) return lp;

  // now, branch on branch_var
  // TODO: could be slightly more efficient by using the current TSP_prob for
  // one of the branches, instead of making two new copies

  // x=0 branch (TODO: should we start with x=1 branch?)
  TSP_prob tp0(*this); // TODO: figure out this constructor bs
  tp0.depth++;
  tp0.fix_var(branch_var, 0.0);
  bool feas0 = tp0.recsolve();
  double objval0 = glp_get_obj_val(tp0.lp) ? feas0 : INF;
  if(feas0) this->ub = min(this->ub, tp0.ub);

  // x=1 branch
  TSP_prob tp1(*this); // TODO: figure out this constructor bs
  tp1.depth++;
  tp1.fix_var(branch_var, 1.0);
  bool feas1 = tp1.recsolve();
  double objval1 = glp_get_obj_val(tp1.lp) ? feas1 : INF;

  if(objval0 < objval1){
    *this = std::move(tp0);
  }else if(feas1){
    *this = std::move(tp1);
  }else{
    return false;
  }

  return true;
}

bool TSP_prob::cp_solve(){
  // solve and add subtour inequalities until none are violated
  int debug = 0;
  while(true){
    printf("Running simplex\n");
    printf("Num rows: %d\n", glp_get_num_rows(this->lp));
    printf("Num cols: %d\n", glp_get_num_cols(this->lp));
    glp_smcp parm;
    glp_init_smcp(&parm);
    //parm.msg_lev = GLP_MSG_OFF;
    assert(glp_simplex(this->lp, &parm) == 0);
    int status = glp_get_status(this->lp);
    if(status != GLP_OPT){
      assert(status == GLP_NOFEAS);
      printf("NOFEAS\n");
      return false;
    }

    if(debug >= 20){
      glp_print_sol(lp, "sol.txt");
      assert(false);
    }

    auto [pmap, weight] = this->min_cut();
    printf("weight is %f\n", weight);
    if(weight >= 2-EPS) break;
    this->add_subtour_constraint(pmap);
    debug++;
  }

  // TODO: look for more cutting planes, e.g. comb inequalities

  return true;
}


pair<TSP_prob::parity_map, double> TSP_prob::min_cut(){
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property<boost::edge_weight_t, int>> undirected_graph;
  typedef boost::graph_traits<undirected_graph>::vertex_descriptor vertex_descriptor;
  typedef boost::property_map<undirected_graph, boost::edge_weight_t>::type weight_map_type;
  typedef boost::property_traits<weight_map_type>::value_type weight_type;

  int ne = this->num_nodes*(this->num_nodes-1)/2;
  /*static*/auto wbuf = new double[ne];

  for(int i = 0; i < ne; i++){
    wbuf[i] = glp_get_col_prim(this->lp, i+1);
  }

  undirected_graph g(edges, edges + ne, wbuf, this->num_nodes, ne);
  BOOST_AUTO(parities, boost::make_one_bit_color_map(num_vertices(g), get(boost::vertex_index, g)));
  double w = boost::stoer_wagner_min_cut(g, get(boost::edge_weight, g), boost::parity_map(parities));

  return {parities, w};
}

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
