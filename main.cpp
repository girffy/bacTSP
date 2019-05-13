#include <stdio.h>
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

vector<point> random_points(int n){
  vector<point> points(n);
  for(int i = 0; i < n; i++){
    points[i] = make_pair(rand01()*100, rand01()*100);
  }

  return points;
}


int main (int argc, char *argv[]){
  int n = 10;
  if(argc > 1){
    n = atoi(argv[1]);
  }
  assert(n > 2);

  vector<point> points = random_points(n);
  //vector<point> points = {{0,0}, {0,1}, {1,0},  {10,10}, {10,11}, {11,10}};
  dmtx d = make_dmtx(points);
  TSP_prob tp(d);
  tp.solve();
  tp.print_wmat();
  tour t = tp.get_tour();
  printf("Optimal tour:\n\n");
  for(int i = 0; i < (int)t.size(); i++){
    printf("%d\t(%f,%f)\n", t[i], points[t[i]].first, points[t[i]].second);
  }
}
