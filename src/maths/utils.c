#include "utils.h"
#include "linalg.h"

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

/*
 * Loads a regressor file where the regressors are in the columns
 * and the samples in the rows. Regressor values should be separated
 * by spaces or tabs.
 * Returns a 2D matrix in the following form: double[regressor][time]
 */
double **rsLoadRegressors(char *path, long *nRegressors, long *nValues, double constantFactor)
{
    FILE *f = fopen(path, "r");
    
    if (f == NULL) {
        fprintf(stderr, "Error: Regressors could not be read.\n");
        return NULL;
    }
    
    /* Read regressor file to get the number of regressors and samples */
    char *line = malloc(sizeof(char)*10000);
    int length = 0;
    *nRegressors = -1L;
    *nValues = 0L;
    
    rewind(f);
    
    long n = 0L;
    double *regressors = NULL;
    while( rsReadline(f, line, &length) ) {
        if ( length < 1 ) {
            continue;
        }
        
        regressors = rsParseRegressorLine(line, &n);
        if (*nRegressors < 0) {
            *nRegressors = n;
        }
        
        if (regressors != NULL) {
            free(regressors);
        }
        *nValues = *nValues+1L;
    }
    
    /* Initialize result matrix */
    rewind(f);
    double **result = d2matrix(*nRegressors, *nValues-1);
    
    /* Fill with constant regressor */
    for ( long v=0L; v<n; v = v+1L ) {
        result[0][v] = constantFactor;
    }
    
    /* Save regressors in result matrix */
    long v=0L;
    while( rsReadline(f, line, &length) ) {
        if ( length < 1 ) {
            continue;
        }
        
        regressors = rsParseRegressorLine(line, &n);
        
        for( long l=0L; l<n; l=l+1L ) {
            result[l+1L][v] = regressors[l];
        }
        
        if (regressors != NULL) {
            free(regressors);
        }
        v=v+1L;
    };
    
    free(line);
    fclose(f);
    
    return result;
}


/*
 * Reads in a single line from a file and returns it.
 */
BOOL rsReadline(FILE *f, char *line, int *length) {
    *length = 0;
    int c;
    
    while(TRUE) {
        c = fgetc(f);
        
        if (c == '\n' || c == '\r' || c == EOF) {
            break;
        }
        
        line[*length] = (char)c;
        *length = *length+1;
    }
    
    line[*length] = '\0';
    
    return c!=EOF;
}

/*
 * Reads in a line from the regressor file and returns the
 * tab- or space-separated values as an ar array of doubles.
 */
double *rsParseRegressorLine(char *line, long *nRegressors) {
    double *regressors;
    char delimiter[] = " \t";
    char *ptr;
    *nRegressors = 0L;
    BOOL endsWithSeparator = TRUE;
    
    size_t lineLength = strlen(line);
    
    // if it doesn't end with a tab or space we need to add one for strtok to work
    if (line[lineLength-1] != '\t' || line[lineLength-1] != ' ') {
        lineLength = lineLength + 1;
        endsWithSeparator = FALSE;
    }
    
    char lineCpy[lineLength+1];
    strcpy(lineCpy, line);
    
    if ( !endsWithSeparator ) {
        lineCpy[lineLength-1] = '\t';
		lineCpy[lineLength]   = '\0';
    }
    
    ptr = strtok(lineCpy, delimiter);
    while (ptr != NULL) {
        ptr = strtok(NULL, delimiter);
        *nRegressors = *nRegressors+1;
    }
    
    regressors = malloc(sizeof(double)*(*nRegressors));
    
    strcpy(lineCpy, line);
    
    if ( !endsWithSeparator ) {
        lineCpy[lineLength-1] = '\t';
    }
    
    long n=0L;
    ptr = strtok(lineCpy, delimiter);
    while(TRUE) {
        if (ptr == NULL) {
            break;
        }
        
        regressors[n] = atof(ptr);
        n = n+1L;
        ptr = strtok(NULL, delimiter);
    }
    
    return regressors;
}