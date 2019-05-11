#include <stdio.h>
#include "TSP_prob.hpp"

double rand01(){
  return abs((float)rand() / (float)RAND_MAX);
}

double l2dist(point p, point q){
  point diff = {p.first-q.first, p.second-q.second};
  return sqrt((diff.first*diff.first + diff.second*diff.second));
}

pair<vector<point>, dmtx> generate_random_TSP(int n){
  vector<point> points(n);
  for(int i = 0; i < n; i++){
    points[i] = make_pair(rand01()*100, rand01()*100);
  }

  vector<vector<double>> d(n);
  for(int i = 0; i < n; i++){
    d[i] = vector<double>(n);
    for(int j = 0; j < n; j++){
      d[i][j] = l2dist(points[i], points[j]);
    }
  }

  return {points, d};
}

int main(){
  int n = 10;
  auto [points, d] = generate_random_TSP(n);
  TSP_prob tp(d);
  tp.solve();
}
