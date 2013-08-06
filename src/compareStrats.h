#ifndef _hyperstrat_RCPP_COMPARESTRATS_H
#define _hyperstrat_RCPP_COMPARESTRATS_H

#include <Rcpp.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */

double rhoij(Rcpp::NumericMatrix qi, Rcpp::NumericMatrix qj, 
	     int k, int l, double n);

double rhoi(Rcpp::NumericMatrix qi, int k, double n);

RcppExport SEXP compareStrats(SEXP q1, SEXP q2);



#endif
