//
//  rsmathutils_internal.h
//  rstools
//
//  Created by André Hoffmann on 5/25/13.
//
//

#include <RInside.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>

#ifndef rstools_rsmathutils_internal_h
#define rstools_rsmathutils_internal_h

Rcpp::NumericVector rsMathCreateVector(double *v, int length);
double* rsMathCovertDoubleVector(Rcpp::NumericVector V);
void rsMathCovertVector(Rcpp::NumericVector V, double *target);

/*Rcpp::NumericMatrix createMatrix(const int n) {
    Rcpp::NumericMatrix M(n,n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            M(i,j) = i*10 + j;
        }
    }
    return(M);
}*/

#endif
