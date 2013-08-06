#include "compareStrats.h"
#include "math.h"

double rhoij(Rcpp::NumericMatrix qi, 
	     Rcpp::NumericMatrix qj, 
	     int k, int l, double n) {
  using namespace Rcpp ;
  double res0 = 0.0;
  for(int s = 0; s < n; ++s) {
    res0 += qi(s,k) * qj(s,l);
  }
  return((1/n) * res0);
}


double rhoi(Rcpp::NumericMatrix qi,  
	     int k, double n) {
  using namespace Rcpp ;
  double res0 = 0.0;
  for(int s = 0; s < n; ++s) {
    res0 += qi(s,k);
  }
  return((1/n) * res0);
}


SEXP compareStrats(SEXP q1, SEXP q2) {
    using namespace Rcpp ;

    NumericMatrix qi(q1); // group assignments
    NumericMatrix qj(q2); // group assignments
    int ni  = qi.ncol();   // num groups for q1
    int nj  = qj.ncol();   // num groups for q2
    double m  = (double)qj.nrow(); // nodes

    double top = 0.0;
    double bottom = 0.0;
    
    for (int k = 0; k < ni; ++k) {
      for (int l = 0; l < nj; ++l) {
	if (rhoij(qi,qj,k,l,m) > 0) {
	  top += (rhoij(qi,qj,k,l,m) * 
		  log( rhoij(qi,qj,k,l,m) / (rhoi(qi,k,m) * rhoi(qj,l,m))));
	}
      }
    }

    for (int k = 0; k < ni; ++k) {
      if (rhoi(qi,k,m) > 0) {
	bottom += -1*(rhoi(qi,k,m) * log(rhoi(qi,k,m)));
      }
    }

    for (int l = 0; l < nj; ++l) {
      if (rhoi(qj,l,m) > 0) {
	bottom += -1*(rhoi(qj,l,m) * log(rhoi(qj,l,m)));
      }
    }

    return(List::create((wrap((2*top) / bottom)), top, bottom));
}

