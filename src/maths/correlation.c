#include "correlation.h"

double rsZCorrelation(const double* X, const double* Y, const size_t length)
{
    /* 
     * This implementation is based on:
     * Weissenbacher, Andreas, et al. 
     * "Correlations and anticorrelations in resting-state functional connectivity MRI: a quantitative comparison of preprocessing strategies."
     * Neuroimage 47.4 (2009): 1408-1416.
     */
    
    const long double N = length;
    const long double Nm1 = length - 1;
    const long double epsilon = 0.0000000001;
    
    // compute means
    long double meanX = 0.0;
    long double meanY = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        meanX = meanX + (long double)X[i] / N;
        meanY = meanY + (long double)Y[i] / N;
    }
    
    // compute standard scores(std. dev)
    long double sX = 0.0;
    long double sY = 0.0;

    for (unsigned int i=0; i<length; i=i+1) {
        sX = sX + powl((long double)X[i]-meanX, 2.0) / Nm1;
        sY = sY + powl((long double)Y[i]-meanY, 2.0) / Nm1;
    }
    
    sX = sqrtl(sX);
    sY = sqrtl(sY);
    
    // compute correlation coeeficient
    long double r = 0.0;
    const long double norm = fmaxl(sX, epsilon) * fmaxl(sY, epsilon) * Nm1;
    
    for (unsigned int i=0; i<length; i=i+1) {
        r = r + (long double)((long double)X[i]-meanX) * ((long double)Y[i]-meanY) / norm;
    }
    
    // Fisher's r-to-z transformation
    const long double half     = 0.5;
    const long double one      = 1.0;
    const long double minusone = -1.0;
    
    r=fminl(r, one-epsilon);
    r=fmaxl(r, minusone+epsilon);
    
    return half * logl( (one + r) / (one - r) );
}

double rsFastZCorrelation(const double* X, const double* Y, const size_t length)
{
    /*
     * This implementation is based on:
     * Weissenbacher, Andreas, et al.
     * "Correlations and anticorrelations in resting-state functional connectivity MRI: a quantitative comparison of preprocessing strategies."
     * Neuroimage 47.4 (2009): 1408-1416.
     */
    
    const double N = length;
    const double Nm1 = length - 1;
    const double epsilon = 0.0000000001;
    
    // compute means
    double meanX = 0.0;
    double meanY = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        meanX = meanX + X[i] / N;
        meanY = meanY + Y[i] / N;
    }
    
    // compute standard scores(std. dev)
    double sX = 0.0;
    double sY = 0.0;
    
    for (unsigned int i=0; i<length; i=i+1) {
        sX = sX + pow(X[i]-meanX, 2.0) / Nm1;
        sY = sY + pow(Y[i]-meanY, 2.0) / Nm1;
    }
    
    sX = sqrt(sX);
    sY = sqrt(sY);
    
    // compute correlation coeeficient
    double r = 0.0;
    const double norm = fmax(sX, epsilon) * fmax(sY, epsilon) * Nm1;
    
    for (unsigned int i=0; i<length; i=i+1) {
        r = r + (X[i]-meanX) * (Y[i]-meanY) / norm;
    }
    
    // Fisher's r-to-z transformation
    const double half     = 0.5;
    const double one      = 1.0;
    const double minusone = -1.0;
    
    r=fmin(r, one-epsilon);
    r=fmax(r, minusone+epsilon);
    
    return half * log( (one + r) / (one - r) );
}

double rsTCorrelation(const double* X, const double* Y, const size_t length)
{
    const double r = rsCorrelation(X, Y, length);
    return r * sqrt( (length-2.0) / (1.0-r*r));
}

double rsCorrelation(const double* X, const double* Y, const size_t length)
{
    return gsl_stats_correlation(X, 1,
                                 Y, 1,
                                 length);
}

double rsMonteCarloZCorrelation(const double* X, const double* Y, const size_t length, const unsigned int repetitions, const unsigned int samplingSize)
{
	double sampledX[samplingSize];
	double sampledY[samplingSize];
	double correlation = 0.0;
	size_t indices[length];
	
	// create random number generator
	gsl_rng_env_setup();
	gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

	// prepare array with indices that will be randomly drawn from
	for (size_t i = 0; i < length; i++) {
		indices[i] = i;
	}

	// repeat drawing samples and computing correlation
	for (unsigned int run=0; run<repetitions; run=run+1) {
		gsl_ran_shuffle(r, indices, length, sizeof(size_t));
		
		//draw samples
		for (unsigned int i=0; i<samplingSize; i=i+1) {
			sampledX[i] = X[indices[i]];
			sampledY[i] = Y[indices[i]];
		}
		
		correlation = correlation + rsZCorrelation(sampledX, sampledY, samplingSize);
	}
	
	gsl_rng_free (r);
	
	return correlation / repetitions;
}