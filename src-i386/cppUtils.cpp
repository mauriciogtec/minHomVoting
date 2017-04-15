//Includes/namespaces
#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

//' @title incidence matrix from adjacency matrix in full form
//' @description takes an adjacency matrix and generates the incidencematrix
//' @param adj an adjacency matrix
//' @return an incidence matrix
//' @export
// [[Rcpp::export]]
IntegerMatrix adj_to_inc(IntegerMatrix adj) {
  int m = adj.nrow();
  CharacterVector v_names(m);
  for (int i=0; i<m; i++) {
    v_names[i] = i + 1;
  }
  int n = m * m;
  CharacterVector e_names(n);
  IntegerMatrix inc(m, n);
  int e_count = 0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      stringstream e;
      e << i + 1 << "_" << j + 1;
      e_names[e_count] = e.str();
      if (adj(i, j) == 1) {
        inc(i, e_count) = 1;
        inc(j, e_count) = -1;
      }
      e_count++;
    }
  }
  inc.attr("dimnames") = List::create(v_names, e_names);
  return inc;
}


//' @title incidence matrix from adjacency matrix in reduced form
//' @description takes an adjacency matrix and generates the incidencematrix
//' @param adj an adjacency matrix
//' @return an incidence matrix
//' @export
// [[Rcpp::export]]
IntegerMatrix adj_to_inc_reduced(IntegerMatrix adj) {
  int m = adj.nrow();
  CharacterVector v_names(m);
  for (int i=0; i<m; i++) {
    v_names[i] = i + 1;
  }
  int n = m * (m - 1) / 2;
  CharacterVector e_names(n);
  IntegerMatrix inc(m, n);
  int e_count = 0;
  for (int i = 0; i < m - 1; i++) {
    for (int j = i + 1; j < m; j++) {
      stringstream e;
      e << i + 1 << "_" << j + 1;
      e_names[e_count] = e.str();
      if (adj(i, j) == 1) {
        inc(i, e_count) = 1;
        inc(j, e_count) = -1;
      }
      e_count++;
    }
  }
  inc.attr("dimnames") = List::create(v_names, e_names);
  return inc;
}