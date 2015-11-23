#include "utils.h"
#include "linalg.h"
#include <sys/types.h>

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
double **rsLoadRegressors(const char *path, long *nRegressors, long *nValues, double constantFactor)
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

double *rsReadRegressorFromStream(FILE *stream, unsigned int *nValues)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    double value;
    unsigned int nBuffer = 10;
    *nValues = 0;
    double *regressor = (double*)rsMalloc(nBuffer * sizeof(double));
    
    while ((read = getline(&line, &len, stream)) != -1) {
        value = atof(line);
        regressor[*nValues] = value;
        *nValues = *nValues + 1;
        
        // Check if we're running out of memory and extend the array if necessary
        if ( *nValues + 1 >= nBuffer ) {
            nBuffer = nBuffer + 10;
            double* tmpRegressor = realloc(regressor, nBuffer * sizeof(double));
            if (tmpRegressor) {
                regressor = tmpRegressor;
            } else {
                fprintf(stderr, "Could not allocate enough memory to save the regressor.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    if (line) free(line);
    
    return regressor;
}

/*
 * Loads a 2D matrix file
 */
double **rsLoadMatrixFromTxt(const char *path, long *nColumns, long *nRows)
{
    FILE *f = fopen(path, "r");

    if (f == NULL) {
        fprintf(stderr, "Error: Matrix could not be read from %s.\n", path);
        return NULL;
    }

    /* Read file to get the number of columns and rows */
    char *line = malloc(sizeof(char)*10000);
    int length = 0;
    *nColumns = 0L;
    *nRows = 0L;

    rewind(f);

    double *columns = NULL;
    while( rsReadline(f, line, &length) ) {
        if ( length < 1 ) {
            continue;
        }

        columns = rsParseRegressorLine(line, nColumns);

        if (columns != NULL) {
            free(columns);
        }
        *nRows = *nRows+1L;
    }

    /* Initialize result matrix */
    rewind(f);
    double **result = d2matrix(*nColumns-1, *nRows-1);

    /* Save columns in result matrix */
    long v=0L;
    while( rsReadline(f, line, &length) ) {
        if ( length < 1 ) {
            continue;
        }

        columns = rsParseRegressorLine(line, nColumns);

        for( long l=0L; l<*nColumns; l=l+1L ) {
            result[l][v] = columns[l];
        }

        if (columns != NULL) {
            free(columns);
        }
        v=v+1L;
    };

    free(line);
    fclose(f);

    return result;
}
