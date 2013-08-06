#include "optimizeColOrder.h"



double getSums(Rcpp::NumericMatrix t1, int a, 
	       Rcpp::NumericMatrix t2, int b) {
    double sum = 0;
    int n = t1.nrow();
    for (int i = 0; i < n; ++i) {
	sum += t1(i,a) * t2(i,b);
    }
    return(sum);
}


SEXP optimizeColOrder(SEXP tab1, SEXP tab2){
    using namespace Rcpp ;
  
    NumericMatrix t1(tab1); // table to align *TO*
    NumericMatrix t2(tab2); // table altering columns
    int ncols = t1.ncol();  // number of columns
    
    NumericVector c1(ncols); // will hold a record of what
    NumericVector c2(ncols); // columns are left as we make selections

    NumericVector res1(ncols); // The pairs of columns are represented 
    NumericVector res2(ncols); // as (res1[i], res2[i])

    // our pairs of columns //
    List l = List(ncols);

    // init our two column vectors //
    for (int i = 0; i < ncols; ++i) {
      c1[i] = i;
      c2[i] = i;
    }

    // traverse all positions of the column vector
    for (int i = 0; i < ncols; ++i) {

	// find the best pairing of whats left
        int n = c1.size(); // what's left
	int nsq = n*n; 
	int m = 0;
	NumericVector colsIndex1(nsq); // keep track of what columns
	NumericVector colsIndex2(nsq); // where used in the sum
	NumericVector as(nsq); // keep track of what columns
	NumericVector bs(nsq); // where used in the sum
	NumericVector sums(nsq);

	if (n > 1) { // if we still have some combinations to check.
	    
	    // first get all the products between the remaining columns
	    for (int a = 0; a < n; ++a) {
		for (int b = 0; b < n; ++b) {
		    colsIndex1[m] = c1[a];  // The actual table indices
		    colsIndex2[m] = c2[b];  // from what's left
		    as[m] = a;              // The index into c1 and c2
		    bs[m] = b;
		    sums[m] = getSums(t1, c1[a], t2, c2[b]);
		    ++m;
		}
	    }

	    // then we need to find the maximum sum ... 
	    // those are the columns to keep
	    double max = 0;
	    int k = 0; // index of best pairing //
	    for (int j = 0; j < nsq; ++j) {
		if (max < sums[j]) {
		    max = sums[j];
		    k = j;
		}
	    }

	    // we save the index of the best pair and 
	    // remove the column index from c1 and c2
	    res1[i] = colsIndex1[k] + 1;
	    res2[i] = colsIndex2[k] + 1; // back to R index land
	    // which c1 & c2 index to erase?
	    c1.erase(as[k]);
	    c2.erase(bs[k]);      
	} else {
	    res1[i] = c1[0] + 1; // back to R index land.
	    res2[i] = c2[0] + 1;	    
	}
    }
    return(wrap(List::create(res1, res2)));
}


