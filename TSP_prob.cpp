#include "TSP_prob.hpp"
#include <assert.h>
#include <stdio.h>

TSP_prob::TSP_prob(dmtx d) : d(d) {
  int n = d.size();
  this->num_nodes = n;
  this->edges = new edge[n*(n-1)/2];
  this->wbuf = new double[n*(n-1)/2];
  this->init_LP();
  this->ia = new int[1 + n*n];
  this->ja = new int[1 + n*n];
  this->ar = new double[1 + n*n];
}

TSP_prob::~TSP_prob(){
  delete this->edges;
  delete this->wbuf;
  delete this->ia;
  delete this->ja;
  delete this->ar;
  glp_delete_prob(this->lp);
}

TSP_prob::TSP_prob(const TSP_prob &other){
  this->lp = glp_create_prob();
  glp_copy_prob(this->lp, other.lp, GLP_ON);

  int n = other.num_nodes;
  int ne = n*(n-1)/2;
  this->num_nodes = n;
  this->d = other.d;
  this->wbuf = new double[ne];
  this->edges = new edge[ne];
  memcpy(this->edges, other.edges, ne * sizeof(edge));

  this->ia = new int[1 + n*n];
  this->ja = new int[1 + n*n];
  this->ar = new double[1 + n*n];
}

TSP_prob& TSP_prob::operator=(TSP_prob&& other){
  this->num_nodes = other.num_nodes;
  this->lp = other.lp;
  this->d = other.d;
  this->edges = other.edges;
  this->wbuf = other.wbuf;
  this->ia = other.ia;
  this->ja = other.ja;
  this->ar = other.ar;

  other.lp = NULL;
  other.edges = NULL;
  other.wbuf = NULL;
  other.ia = NULL;
  other.ja = NULL;
  other.ar = NULL;

  return *this;
}

void TSP_prob::init_LP(){
  int n = this->d.size();
  glp_prob *lp;
  lp = glp_create_prob();
  glp_set_prob_name(lp, "TSP");
  glp_set_obj_dir(lp, GLP_MIN);
  glp_add_rows(lp, n);
  glp_add_cols(lp, n*(n-1)/2);

  // add rows/constraints
  char buffer[50];
  for(int i=0; i<n; i++){
    sprintf(buffer, "degcon%d", i);
    glp_set_row_name(lp, i+1, buffer);
    glp_set_row_bnds(lp, i+1, GLP_FX, 2.0, 2.0);
  }

  // add variables
  int vidx = 1;
  for(int i=0; i<n; i++){
    for(int j=i+1; j<n; j++){
      sprintf(buffer, "x%d,%d", i, j);
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
      this->ia[aidx] = a+1;
      this->ja[aidx] = vidx;
      this->ar[aidx] = 1.0;
      aidx++;

      this->ia[aidx] = b+1;
      this->ja[aidx] = vidx;
      this->ar[aidx] = 1.0;
      aidx++;
      vidx++;
    }
  }

  glp_load_matrix(lp, aidx-1, this->ia, this->ja, this->ar);
}

void TSP_prob::add_subtour_constraint(parity_map pmap){
  int n = this->num_nodes;
  // TODO: make ind and val scratch variables in TSP_prob
  int *ind = new int[1 + n];
  double *val = new double[1 + n];
  int eidx = 1;
  glp_add_rows(this->lp, 1);
  int r = glp_get_num_rows(this->lp);
  glp_set_row_name(this->lp, r, "subtour");
  glp_set_row_bnds(this->lp, r, GLP_LO, 2.0, 0.0);

  for(int src = 0; src < n; src++){
    for(int dest = src+1; dest < n; dest++){
      if(get(pmap, src) != get(pmap, dest)){
        edge e = {src, dest};
        int idx = edge_to_var(e);
        ind[eidx] = idx;
        val[eidx] = 1.0;
        eidx++;
      }
    }
  }

  glp_set_mat_row(this->lp, r, eidx-1, ind, val);

  delete ind;
  delete val;
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

int TSP_prob::edge_to_var(int src, int dest){
  assert(src < dest);
  int naive = this->num_nodes*src - dest;
  return naive - src*(src+1)/2 + 1; // +1 for glpk
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
  while(true){
    assert(glp_simplex(this->lp, NULL) == 0);
    int status = glp_get_prim_stat(this->lp);
    if(status == GLP_INFEAS) return false;

    auto [pmap, weight] = this->min_cut();
    if(weight >= 2-EPS) break;
    this->add_subtour_constraint(pmap);
  }

  // TODO: look for more cutting planes, e.g. comb inequalities

  return true;
}


pair<TSP_prob::parity_map, double> TSP_prob::min_cut(){
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property<boost::edge_weight_t, int> > undirected_graph;
  typedef boost::graph_traits<undirected_graph>::vertex_descriptor vertex_descriptor;
  typedef boost::property_map<undirected_graph, boost::edge_weight_t>::type weight_map_type;
  typedef boost::property_traits<weight_map_type>::value_type weight_type;

  int ne = this->num_nodes*(this->num_nodes-1)/2;

  for(int i = 0; i < ne; i++){
    wbuf[i] = glp_get_col_prim(this->lp, i+1);
  }

  undirected_graph g(edges, edges + ne, wbuf, this->num_nodes, ne);
  BOOST_AUTO(parities, boost::make_one_bit_color_map(num_vertices(g), get(boost::vertex_index, g)));
  double w = boost::stoer_wagner_min_cut(g, get(boost::edge_weight, g), boost::parity_map(parities));

  return {parities, w};
}
