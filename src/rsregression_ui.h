#ifndef rstools_rsregression_ui_h
#define rstools_rsregression_ui_h

#include <stdio.h>
#include <strings.h>
#include <regex.h>
#include <glib.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"
#include "src/utils/rsui.h"

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

    long nRegressors;
    long nRegressorValues;

    double **regressors;
    double ***mask;

    rsUIInterface *interface;
    
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
void rsRegressionBuildInterface(rsRegressionParameters *p);
void rsRegressionFreeParams(rsRegressionParameters* p);

#ifdef __cplusplus
}
#endif

#endif
