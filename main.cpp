#include <stdio.h>
#include <fstream>
#include "TSP_prob.hpp"

double rand01(){
  return abs((float)rand() / (float)RAND_MAX);
}

double l2dist(point p, point q){
  point diff = {p.first-q.first, p.second-q.second};
  return sqrt((diff.first*diff.first + diff.second*diff.second));
}

dmtx make_dmtx(vector<point> points){
  int n = points.size();
  vector<vector<double>> d(n);
  for(int i = 0; i < n; i++){
    d[i] = vector<double>(n);
    for(int j = 0; j < n; j++){
      d[i][j] = l2dist(points[i], points[j]);
    }
  }

  return d;
}

// expects the first line to be n, the number of cities, followed by a matrix of
// n^2 entries indicating distances
dmtx read_dmtx(char *fname){
  ifstream inp;
  inp.open(fname);
  int n;
  inp >> n;
  dmtx d(n);
  for(int i = 0; i < n; i++){
    d[i] = vector<double>(n);
    for(int j = 0; j < n; j++){
      assert(inp >> d[i][j]);
    }
  }

  assert(!(inp >> n));

  return d;
}

vector<point> random_points(int n){
  vector<point> points(n);
  for(int i = 0; i < n; i++){
    points[i] = make_pair(rand01()*100, rand01()*100);
  }

  return points;
}


int main (int argc, char *argv[]){
  int n;
  vector<point> points;
  dmtx d;
  if(argc <= 1){
    printf("Usage: %s <num_cities> | %s file <dmtx>\n", argv[0], argv[0]);
    exit(1);
  }else if(strcmp(argv[1], "file") == 0){
    assert(argc == 3);
    d = read_dmtx(argv[2]);
  }else{
    assert(argc == 2);
    n = atoi(argv[1]);
    assert(n > 2);
    vector<point> points = random_points(n);
    d = make_dmtx(points);
  }

  TSP_prob tp(d);
  tp.solve();
  tp.print_wmat();
  tour t = tp.get_tour();
  printf("Optimal tour:\n\n");
  for(int i = 0; i < (int)t.size(); i++){
    if(points.size() == 0){
      printf("%d\n", t[i]);
    }else{
      printf("%d\t(%f,%f)\n", t[i], points[t[i]].first, points[t[i]].second);
    }
  }
  printf("tour length: %f\n", glp_get_obj_val(tp.lp));
}
