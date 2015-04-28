#include <Rcpp.h>
using namespace Rcpp;

// really vectorized power function (exponent)
NumericVector vec_pow(NumericVector base, NumericVector exp) {
  NumericVector out = clone(base);
  int n = out.length();
  for(int i = 0; i < n; i++){
    out[i] = pow(base[i], exp[i]);
  }
  return out;
}

// simpler range construction
NumericVector range(int start, int end){
  NumericVector x;
  for( int i = start; i <= end; i++ ) x.push_back( i );
  return x;    
}

// logistic function
double lf(double x) {
  return exp(x)/(1+exp(x));
}

// Three Paramater Logistic (MultiDim) (Segall, 1996)
// [[Rcpp::export]]
List PROB_3PLM(NumericVector theta, NumericMatrix a, NumericVector b, NumericVector c, NumericVector u) {
  int Q = a.ncol(), K = a.nrow(), M = 2;
  bool posterior = is_false(any(is_na(u)));
  NumericVector aux(K), l(K), d(K), D(K);
  NumericMatrix P(K, M);
  
  // loop over items
  for (int k = 0; k < K; k++){
    aux[k] = 0;
    
    // loop over dimensions
    for (int q = 0; q < Q; q++){
      aux[k] += -a(k,q) * (theta[q] - b[k]);
    }
    
    P(k,1) = c[k] + (1-c[k])/(1+exp(aux[k]));
    P(k,0) = 1 - P(k,1);
  }
  
  if(posterior){
    NumericVector p = P(_,1), p2 = pow(p,2), q = P(_,0);
    
    
    l = vec_pow(p, u) * vec_pow(q, 1 - u);
    d = ((p-c) * (u-p)) / ((1-c) * p);
    D = (q *(p-c)*(c*u-p2)) / (p2*pow(1-c,2));
  }
  
  List out = List::create(_["P"] = P,
  _["l"] = l,
  _["d"] = d,
  _["D"] = D);
  
  return out;
}

// graded response model (Samejima, 1969)
// [[Rcpp::export]]
List PROB_GRM(NumericVector theta, NumericMatrix a, NumericMatrix b, NumericVector u) {
  int K = a.nrow(), M = b.ncol() + 1, m = 0;
  double at = 0;
  bool posterior = is_false(any(is_na(u)));
  
  NumericVector l(K), d(K), D(K);
  NumericMatrix P(K, M);
  for (int i = 0; i < P.ncol() * P.nrow(); i++) P[i] = NA_REAL;
  
  // loop over items
  for(int k = 0; k < K; k++){
    // no of response categories
    m = sum(!is_na(b(k,_)));
    
    // inner product alpha theta
    at = sum((a(k,_) * theta));
    
    NumericVector Psi(m+2, 1.0); // default to 1 (first element)
    // loop over response categories.
    for (int i = 1; i < m + 2; i++){
      if (i <= m) {
        Psi[i] = lf((at - b(k, i-1)));
      } else { 
        // last element is 0
        Psi[i] = 0;
      }
      
      // build probability matrix.
      P(k,i-1) = Psi[i - 1] - Psi[i];
    }
    
    // build likelihood and derivatives
    if(posterior){
      l[k] = P(k,u[k]);
      d[k] = 1 - Psi[u[k]] - Psi[u[k]+1];
      D[k] = -(Psi[u[k]] * (1-Psi[u[k]]) + Psi[u[k]+1] * (1-Psi[u[k]+1]));
    }
  }
  
  List out = List::create(_["P"] = P,
  _["l"] = l,
  _["d"] = d,
  _["D"] = D);
  
  return out;
}

//  Sequential Model (Tutz, 1990)
// [[Rcpp::export]]
List PROB_SM(NumericVector theta, NumericMatrix a, NumericMatrix b, NumericVector u) {
  int K = a.nrow(), M = b.ncol() + 1, m = 0;
  double at = 0;
  bool posterior = is_false(any(is_na(u)));
  
  NumericVector l(K), d(K), D(K);
  NumericMatrix P(K, M);
  for (int i = 0; i < P.ncol() * P.nrow(); i++) P[i] = NA_REAL;
  
  // loop over items
  for(int k = 0; k < K; k++){
    // no of response categories
    m = sum(!is_na(b(k,_)));
    
    // inner product alpha theta
    at = sum((a(k,_) * theta));
    
    NumericVector Psi(m+2, 1.0); // default to 1 (first element)
    double aux = 1.0; // keeping track of product so far.
    
    // loop over response categories.
    for (int i = 1; i < m + 2; i++){
      if (i <= m) {
        Psi[i] = lf((at - b(k, i-1)));
      } else { 
        // last element is 0
        Psi[i] = 0;
      }
      
      // build probability matrix.
      aux *= Psi[i-1];
      P(k, i-1) = aux * (1 - Psi[i]);
    }
    
    if (posterior){
      l[k] = P(k,u[k]);
      
      aux = 0;
      D[k] = 0;
      for (int i = 0; i <= u[k]; i++) {
        aux += 1 - Psi[i];
        D[k] -= Psi[i] * (1-Psi[i]);
      }
      
      d[k] = aux - Psi[u[k]+1];
      D[k] -= Psi[u[k]+1] * (1-Psi[u[k]+1]);
    }
  }
  
  List out = List::create(_["P"] = P,
  _["l"] = l,
  _["d"] = d,
  _["D"] = D);
  
  return out;
}


// Generalised Partial Credit Model (Muraki, 1992)
// [[Rcpp::export]]
List PROB_GPCM(NumericVector theta, NumericMatrix a, NumericMatrix b, NumericVector u) {
  int K = a.nrow(), M = b.ncol() + 1, m = 0;
  double at = 0;
  bool posterior = is_false(any(is_na(u)));
  
  NumericVector l(K), d(K), D(K);
  NumericMatrix P(K, M);
  for (int i = 0; i < P.ncol() * P.nrow(); i++) P[i] = NA_REAL;
  
  // loop over items
  for(int k = 0; k < K; k++) {    
    // no of response categories
    m = sum(!is_na(b(k,_)));
    
    // inner product alpha theta
    at = sum((a(k,_) * theta));
    
    // loop over response categories
    NumericVector aux(m);
    for (int i = 0; i < m; i++){
      aux[i] = exp((i+1) * at - b(k,i));
    }
    
    double aux2 = sum(aux), remainder = 1.0;
    for (int i = 0; i < m; i++){
      P(k, i+1) = aux[i] / (1 + aux2);
      remainder -= P(k, i+1);
    }
    P(k, 0) = remainder;
  
    if (posterior){
      NumericVector mi = range(1,m);
      NumericMatrix pi = P(Range(k,k), Range(1,m));
      double mp = sum(mi*pi);
      
      l[k] = P(k,u[k]);
      d[k] = u[k] - mp;
      D[k] = -sum((mi * pi) * (mi - mp));
    }
  }
  
  List out = List::create(_["P"] = P,
  _["l"] = l,
  _["d"] = d,
  _["D"] = D);
  
  return out;
}