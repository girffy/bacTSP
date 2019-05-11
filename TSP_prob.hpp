#include <vector>
#include <utility>
#include <cstdlib>
#include <glpk.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/one_bit_color_map.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/typeof/typeof.hpp>

using namespace std;
typedef pair<double, double> point;
typedef pair<int, int> edge;
typedef vector<vector<double>> dmtx;
typedef vector<int> tour;

class TSP_prob {
  public:
    int num_nodes;
    glp_prob *lp;
    dmtx d;
    edge *edges;
    double *wbuf; // a weight buffer used for min cut

    TSP_prob(dmtx d);
    ~TSP_prob();
    TSP_prob(const TSP_prob &other);
    TSP_prob& operator=(TSP_prob&& other);

  // solve the TSP instance; this->lp will have the solution as its primal
  // variable
  void solve();

  private:
    // a type for the way boost produces a graph cut
    typedef boost::one_bit_color_map<boost::vec_adj_list_vertex_id_map<boost::no_property, long unsigned int> >
        parity_map;

    /* some state variables used during the course of solve() */
    double lb; // the best known lower bound for the problem
    double ub; // the best known upper bound for the problem
    int depth; // the current branching depth

    /* 3 scratch arrays used for interfacing with GLPK */
    int *ia;
    int *ja;
    double *ar;

    /* methods for massaging LPs */
    // initialize the lp object of a TSP_prob; also initialize the edges
    void init_LP();

    // add a subtour constraint for the given cut to the TSP LP
    void add_subtour_constraint(parity_map pmap);

    // add a constraint forcing the given variable to have a specific value
    void fix_var(int var, double value);

    // given an edge (src,dest), determines the index of its column in the LP
    int edge_to_var(edge e);
    int edge_to_var(int src, int dest);

    // given a LP column index i, determines the edge it corresponds to
    edge var_to_edge(int i);



    // recursively solve this TSP_prob via branch-and-cut. The TSP problem is
    // described by this->lp, which may have additional cutting plane and
    // variable-fixing constraints that this function will respect. This
    // function will compute the best tour that satisfies the constraints of
    // this->lp which has length less than ub, or return false if no such tour
    // exists
    bool recsolve();

    // solve the given LP, adding cutting planes as necessary, at least to ensure no
    // subtour constraints are violated, and ideally try as hard as possible to find
    // other cutting planes to refine fractional solutions. Return false if adding
    // valid cutting planes made the problem infeasible
    bool cp_solve();

    // extract a weighted undirected graph from current solution, and find a
    // minimum-weight cut for it. Return a parity map, mapping vertices to
    // true/false according to their side of the cut, and the weight of the cut
    pair<parity_map, double> min_cut();
};
