/*******************************************************************
 *
 * rsmathutils.c
 *
 * 
 * André Hoffmann
 *******************************************************************/

#include "utils.h"

static gsl_rng *rsRandomNumberGenerator = NULL;

double rsSampleSineWave(const double sampling_rate, const double f, const int t)
{
    return sin((double)t * 2.0 * M_PI * f * sampling_rate);
}

double rsSampleCosineWave(const double sampling_rate, const double f, const int t)
{
    return cos((double)t * 2.0 * M_PI * f * sampling_rate);
}

double rsSigmoidRolloff(const double nBins, const double rolloff, const double bin)
{
    return rsSigmoid(rolloff/nBins, bin);
}

double rsSigmoid(const double rolloff, const double x)
{
    return 1.0 - tanh(rolloff*M_PI*pow(x, 2.0));
}

int rsCountDigits(int n)
{
    return (n == 0) ? 1 : floor(log10(abs(n))) + 1;
}

gsl_rng *rsGetRandomNumberGenerator()
{
	if ( ! rsRandomNumberGenerator ) {
		gsl_rng_env_setup();
		rsRandomNumberGenerator = gsl_rng_alloc(gsl_rng_default);
	}
	
	return rsRandomNumberGenerator;
}

void rsDestroyRandomNumberGenerator()
{
	gsl_rng_free(rsRandomNumberGenerator);
	rsRandomNumberGenerator = NULL;
}