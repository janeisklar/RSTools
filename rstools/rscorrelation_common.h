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
#include "rsniftiutils.h"
#include "rsmathutils.h"

struct rsCorrelationParameters {
    char *inputpath;
	char *maskpath;
	char *outputpath;
    char *savemaskpath;
    char *saveBetasPath;
    char *saveResidualsPath;
    char *saveFittedPath;
	
	short xDim;
    short yDim;
    short zDim;
    short vDim;
    
	short pixtype;
	size_t dt;
    
    float inter;
    float slope;
    
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