//
//  rsregression_common.h
//  rstools
//
//  Created by André Hoffmann on 6/27/13.
//
//

#ifndef rstools_rsregression_common_h
#define rstools_rsregression_common_h

#include <stdio.h>
#include <strings.h>
#include <regex.h>
#include <nifti1.h>
#include <fslio.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"

struct rsRegressionParameters {
    char *inputpath;
	char *maskpath;
	char *regressorspath;
	char *regressorCommentPath;
	char *regressorComment;
    char *savemaskpath;
    char *saveBetasPath;
    char *saveResidualsPath;
    char *saveFittedPath;
	
	short xDim;
    short yDim;
    short zDim;
    short vDim;
    
	size_t pixtype;
	short dt;
    
    float inter;
    float slope;
    
    double f1;
    double f2;
    double sampling_rate;
    
    BOOL verbose;
    BOOL filterActive;
    BOOL parametersValid;
    BOOL zScoreRegression;
    
    long nRegressors;
    long nRegressorValues;
    
    double **regressors;
    
    FSLIO *fslio;
    double ***mask;
    
    double nyquist_frequency;
    double bin_width;
    
    int nFrequencyBinsLow;
    int nFrequencyBinsHigh;
    int nFrequencyBins;
    int nFrequencyRegressors;
    
    double *frequencyBins;
    int nAllRegressors;
    double **allRegressors;
    
    int threads;
    size_t wordsize;
};

BOOL rsReadline(FILE *f, char *line, int *length);
double **rsLoadRegressors(char *path, long *nRegressors, long *nValues, double constantFactor);
double *rsParseRegressorLine(char *line, long *nRegressors);
void rsRegressionPrintHelp();
struct rsRegressionParameters rsRegressionLoadParams(int argc, char * argv[]);

#endif
