#ifndef rstools_rsregression_ui_h
#define rstools_rsregression_ui_h

#include <stdio.h>
#include <strings.h>
#include <regex.h>
#include <nifti1.h>
#include <fslio.h>
#include <glib.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *inputpath;
    char *maskpath;
    char *regressorspath;
    char *regressorCommentPath;
    char *regressorComment;
    char *savemaskpath;
    char *saveBetasPath;
    char *saveResidualsPath;
    char *saveFittedPath;
    char *callString;
    char *comment;
    
    double freqLow;
    double freqHigh;
    double TR;
    
    BOOL verbose;
    BOOL filterActive;
    BOOL parametersValid;
    BOOL zScoreRegression;
    
    rsNiftiFile *input;
    rsNiftiFile *betas;
    rsNiftiFile *residuals;
    rsNiftiFile *fitted;

    GOptionContext *context;

    long nRegressors;
    long nRegressorValues;

    double **regressors;
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
    rsReportProgressCallback *progressCallback;

} rsRegressionParameters;

rsRegressionParameters* rsRegressionParseParams(int argc, char * argv[]);
rsRegressionParameters* rsRegressionInitParameters();
void rsRegressionFreeParams(rsRegressionParameters* p);
void rsRegressionPrintHelp(rsRegressionParameters* p);

#ifdef __cplusplus
}
#endif

#endif
