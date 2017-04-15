//Includes/namespaces
#include <Rcpp.h>
#include <iostream>
#include <list>
using namespace Rcpp;
using namespace std;

//' @title find edges
//' @description takes an adjacency matrix and generates a list of edges
//' @return a list with edges
//' @export
// [[Rcpp::export]]
IntegerMatrix find_edges(IntegerMatrix adj) {
  int n = adj.nrow();
  list<int> x;
  list<int> y;
  int n_edges = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      if (adj(i, j) == 1) {
        x.push_back(i + 1);
        y.push_back(j + 1);
        n_edges++;
      }
    }
  }
  IntegerMatrix edges(n_edges, 2);
  for (int i = 0; i < n_edges; i++) {
    edges(i, 0) = x.front();
    edges(i, 1) = y.front();
    x.pop_front();
    y.pop_front();
  }
  return edges;
}

//' @title incidence matrix from adjacency matrix in reduced form
//' @description takes an adjacency matrix and generates the incidencematrix
//' @param adj an adjacency matrix
//' @details the incidence matrix \code{E} has entry \code{e_{ij}} equal to
//' \itemize{
//'     \code{1} if an edge starts in node \code{i} and finishes in \code{j}
//'     \code{-1} if an edge ends in node \code{i} and finishes in \code{j}
//'     \code{0} if nodes \code{i} and \code{j} are not connected
//' }
//' @return an incidence matrix
//' @export
// [[Rcpp::export]]
IntegerMatrix adj_to_inc(IntegerMatrix adj) {
  int m = adj.nrow();
  CharacterVector v_names(m);
  for (int i=0; i<m; i++) {
    v_names[i] = i + 1;
  }
  IntegerMatrix edges = find_edges(adj);
  int n = edges.nrow();
  CharacterVector e_names(n);
  IntegerMatrix inc(m, n);
  int i,j;
  for (int k = 0; k < n; k++) {
    stringstream e;
    i = edges(k, 0) - 1;
    j = edges(k, 1) - 1;
    e << i + 1 << "_" << j + 1;
    e_names[k] = e.str();
    inc(i, k) = 1;
    inc(j, k) = -1;
  }
  inc.attr("dimnames") = List::create(v_names, e_names);
  return inc;
}