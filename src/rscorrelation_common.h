//
//  rscorrelation_common.h
//  rstools
//
//  Created by André Hoffmann on 7/9/13.
//
//

#ifndef rstools_rscorrelation_common_h
#define rstools_rscorrelation_common_h

#include <stdio.h>
#include <nifti1.h>
#include <fslio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"

#define RSTOOLS_CORRELATION_CONVERSION_NONE 1
#define RSTOOLS_CORRELATION_CONVERSION_Z    2
#define RSTOOLS_CORRELATION_CONVERSION_T    3

struct rsCorrelationParameters {
    char *inputpath;
	char *maskpath;
	char *outputpath;
	char *commentpath;
	char *comment;
    char *savemaskpath;
    char *saveBetasPath;
    char *saveResidualsPath;
    char *saveFittedPath;
	
	short xDim;
    short yDim;
    short zDim;
    short vDim;
    
    short delay;
    
	short pixtype;
	size_t dt;
    
    float inter;
    float slope;
    
    short conversionMode;

	unsigned int monteCarloRepetitions;
	unsigned int monteCarloSampleSize;
    
    BOOL verbose;
    BOOL parametersValid;
    
    unsigned int nRegressorValues;
    double *regressor;
    
    FSLIO *fslio;
    FSLIO *fslioCorrelation;
    double ***mask;
    double ***correlation;
    
    int threads;
    size_t wordsize;
};

void                            rsCorrelationPrintHelp();
struct rsBandpassParameters     rsBandpassLoadParams(int argc, char * argv[]);
struct rsCorrelationParameters  rsCorrelationLoadParams(int argc, char * argv[]);
void                            rsCorrelationWriteCorrelationFile(struct rsCorrelationParameters* p);
void                            rsCorrelationFree(struct rsCorrelationParameters* p);
double*                         rsReadRegressorFromStandardInput(unsigned int *nValues);

#endif