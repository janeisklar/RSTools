#include "stats.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// See: http://en.wikipedia.org/wiki/Student's_t-test#One-sample_t-test                                     //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double rsOneSampleTTest(const double *data, const unsigned int length, const double mu)
{
    // compute mean
    double mean = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        mean = mean + data[i] / length;
    }
    
    // compute standard scores(std. dev)
    double std = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        std = std + pow(data[i]-mean, 2.0) / ((double)length-1.0);
    }
    
    std = sqrt(std);
    
    return (mean - mu) / (std / sqrt(length));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following two functions are taken from:                                                              //
// http://www.spraak.org/documentation/doxygen/src/lib/math/erfinv.c/view                                   //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*************************************************************
 *  __  _   _   _   _                                        *
 * (_  |_) |_) |_| |_| |_/  Speech Processing, Recognition & *
 * __) |   | \ | | | | | \  Automatic Annotation Kit         *
 *                                                           *
 * Copyright 2006,2007,2008 K.U.Leuven                       *
 *                                                           *
 * Use of this software is governed by a License Agreement.  *
 * Modifications should be properly credited in line with    *
 * the License Agreement.                                    *
 *************************************************************/

/*
 * Calculate the inverse error function.
 * The fast version is only correct up to 6 digits.
 */
float rsFastErfInv(float x) {
    float tmp;
    int   neg;
    
    if((neg=(x < 0.0)))
        x = -x;
    if(x <= 0.7) {
        tmp = x*x;
        x *= (((-0.140543331*tmp+0.914624893)*tmp-1.645349621)*tmp+0.886226899)/((((0.012229801*tmp-0.329097515)*tmp+1.442710462)*tmp-2.118377725)*tmp+1.0);
    } else {
        tmp = sqrt(-log(0.5*(1.0-x)));
        x = (((1.641345311*tmp+3.429567803)*tmp-1.624906493)*tmp-1.970840454)/((1.637067800*tmp+3.543889200)*tmp+1.0);
    }
    return(neg?-x:x);
}

/*
 * Calculate the inverse error function.
 * Uses fast_erfinv and performs a two steps of Newton-Raphson correction to
 * achieve full accuracy.
 */
double rsErfInv(const double x) {
    const double tmp = rsFastErfInv(x);
    return tmp - (erf(tmp)-x)*exp(tmp*tmp)*0.886226925452757941 - (erf(tmp)-x)*exp(tmp*tmp)*0.886226925452757941;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following two functions are taken from:                                                              //
// https://gist.github.com/timflutre/1784199                                                                //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Resolve a sequence of ties.
 *
 * The input ranks array is expected to take the same value for all indices in
 * tiesTrace. The common value is recoded with the average of the indices. For
 * example, if ranks = <5,8,2,6,2,7,1,2> and tiesTrace = <2,4,7>, the result
 * will be <5,8,3,6,3,7,1,3>.
 *
 * Source: http://commons.apache.org/math/apidocs/src-html/org/apache/commons/math/stat/ranking/NaturalRanking.html#line.312
 */
void rsRankingResolveTies(double *ranks, const size_t *tiesTrace, const size_t n_ties) {
    size_t i;
    
    // constant value of ranks over tiesTrace
    const double c = ranks[tiesTrace[0]];
    
    // new rank (ie. the average of the current indices)
    const double avg = (2*c + n_ties - 1) / 2;
    
    for(i=0; i<n_ties; ++i) {
        ranks[tiesTrace[i]] = avg;
    }
}

/*
 * Rank data using the natural ordering on doubles, ties being resolved by taking their average.
 *
 * Source: http://commons.apache.org/math/apidocs/src-html/org/apache/commons/math/stat/ranking/NaturalRanking.html#line.190
 */
void rsSpearmannRank(double *ranks, const double *data, const size_t n) {
    size_t i;
    double *d = malloc(sizeof(double)*n);
    size_t *p = malloc(sizeof(size_t)*n);
    
    // copy the input data and sort them
    for(i=0; i<n; ++i)
        d[i] = data[i];
    gsl_sort(d, 1, n);
    
    // get the index of the input data as if they were sorted
    gsl_sort_index(p, data, 1, n);
    
    // walk the sorted array, filling output array using sorted positions, resolving ties as we go
    size_t pos = 1;
    ranks[p[0]] = pos;
    size_t n_ties = 1;
    size_t * tiesTrace = (size_t*) calloc (1, sizeof(size_t));
    tiesTrace[0] = p[0];
    for(i=1; i<n; ++i) {
        if(d[i] - d[i-1] > 0) {
            pos = i + 1;
            if(n_ties > 1) {
                rsRankingResolveTies(ranks, tiesTrace, n_ties);
            }
            tiesTrace = (size_t*) realloc(tiesTrace, sizeof(size_t));
            n_ties = 1;
            tiesTrace[0] = p[i];
        } else {
            ++n_ties;
            tiesTrace = (size_t*) realloc(tiesTrace, n_ties * sizeof(size_t));
            tiesTrace[n_ties-1] = p[i];
        }
        ranks[p[i]] = pos;
    }
    if(n_ties > 1) {
        rsRankingResolveTies(ranks, tiesTrace, n_ties);
    }
    
    free(tiesTrace);
    free(d);
    free(p);
}