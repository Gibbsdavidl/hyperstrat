#include "stratem.h"


Rcpp::NumericMatrix init(int n, int ng) {
  using namespace Rcpp ;

  NumericMatrix x = NumericMatrix(n,ng);
  NumericMatrix q = NumericMatrix(n,ng);
  NumericVector z;

  for (int i=0; i < n; ++i) {
    for (int r=0; r < ng; ++r) {
      z = runif(1);
      x(i,r) = z[0];
    }
  }

  for (int i=0; i < n; ++i) {
    for (int r=0; r < ng; ++r) {
      q(i,r) = x(i,r) / sum(x.row(i));
    }
  }
  return(q);
}



Rcpp::NumericMatrix updateTheta(Rcpp::NumericMatrix a, 
				Rcpp::NumericMatrix q,
				int n, int m, int ng) {
  using namespace Rcpp ;
  Rcpp::NumericMatrix theta(ng, m);
  for (int r = 0; r < ng; ++r) {
    for (int j = 0; j < m; ++j) {
      theta(r,j) = Rcpp::sum(q(_,r) * a(j,_)) / Rcpp::sum(q(_,r));
    }
  }
  return(theta);
}


double updateQ(Rcpp::NumericMatrix a, 
	       Rcpp::NumericMatrix theta, 
	       int i, int r, int ng, int m) {
  // So here, what is the probability over all edges
  // for node i being in group r.... compared to
  // the sum over all groups.
  using namespace Rcpp ;
  double top = 1;
  for (int j = 0; j < m; ++j) { // for each edge in group r
    top *= pow(theta(r,j), a(j,i)) * pow((1-theta(r,j)), (1-a(j,i)));
  }
  double bot = 0;
  for (int g = 0; g < ng; ++g) { // all groups
    double x = 1;
    for (int j = 0; j < m; ++j) {
      x *= pow(theta(g,j), a(j,i)) * pow((1-theta(g,j)), (1-a(j,i)));
    }
    bot += x;
  }
  return(top/bot);
}


Rcpp::NumericMatrix updateQs(Rcpp::NumericMatrix a, 
			     Rcpp::NumericMatrix theta, 
			     int n, int ng, int m) {
  using namespace Rcpp ;
  Rcpp::NumericMatrix q(n, ng);
  for (int r = 0; r < ng; ++r) {
    for (int i = 0; i < n; ++i) {
      q(i,r) = updateQ(a,theta,i,r,ng,m);
    }
  }
  return(q);
}


double expectation(Rcpp::NumericMatrix q,
		   Rcpp::NumericMatrix a,
		   Rcpp::NumericMatrix theta,
		   int n, int ng, int m) {
  using namespace Rcpp ;
  double e = 0.0;
  double x = 0.0;
  double y = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int r = 0; r < ng; ++r) {
      for (int j = 0; j < m; ++j) {
	if (theta(r,j) != 0) {
	  x = std::log(theta(r,j));
	} else {
	  x = -1000;
	}
	if ((1-theta(r,j)) != 0) {
	  y = std::log(1-theta(r,j));
	} else {
	  y = -1000;
	}
	// prob of indiv i in group r times
	// the prob of group r having attribute j
	// plus prob of group r not having attribute j
	e = e + q(i,r) * (a(j,i) * x + (1-a(j,i)) * y);
      }
    }
  }
  return(e);
}


SEXP stratem(SEXP dat, SEXP groups, SEXP epsilon, SEXP maxiter){
    using namespace Rcpp ;

    NumericMatrix a(dat); // adjacency matrix
    NumericVector ngroups(groups);
    NumericVector maxi(maxiter); // maximum iterations
    NumericVector eps(epsilon);
    
    int ng = ngroups[0]; // number of groups
    int ni = maxi[0];    // number of likelihood iterations
    int n  = a.ncol();   // nodes
    int m  = a.nrow();   // edges

    int    epoch = 0;
    double e = eps[0];
    double Lold = 0.0;
    double Lnew = 0.0;

    // q is probability of inidiviual i in group r
    NumericMatrix q = init(n, ng); 

    // theta is probability for group r in edge j
    NumericMatrix theta = updateTheta(a, q, n, m, ng);

    // Initial Expectation
    Lnew = expectation(q,a,theta,n,ng,m);

    while((std::abs(Lnew - Lold) > e) && (ni > epoch)) {

      // Maximization
      q = updateQs(a, theta, n, ng, m);
      theta = updateTheta(a, q, n, m, ng);
      
      // Expectation
      Lold = Lnew;
      Lnew = expectation(q,a,theta,n,ng,m);
    
      epoch += 1;
      //printf("epoch: %d\tL: %f\n", epoch, Lnew);
    }
    if (ni > 0) {
      return(wrap(List::create(q,theta,Lnew)));
    } else {
      return(NULL);
    }
}

