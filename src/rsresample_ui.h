#ifndef rstools_resample_ui_h
#define  rstools_resample_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *inputpath;
    char *outputpath;
    char *callString;
    char **regressors;
    char **regressorInputs;
    char **regressorOutputs;
    
    double inputTR;
    double outputTR;
    
    int order;
    int nRegressors;
    
    BOOL verbose;
    BOOL parametersValid;
    
    rsNiftiFile *input;
    rsNiftiFile *output;

    rsUIInterface *interface;

    int threads;
    
    rsReportProgressCallback *progressCallback;
    
} rsResampleParameters;

rsResampleParameters *rsResampleParseParams(int argc, char * argv[]);
rsResampleParameters *rsResampleInitParameters();
void rsResampleFreeParams(rsResampleParameters *p);
void rsResampleBuildInterface(rsResampleParameters *p);
BOOL rsResampleParseRegressorPath(char **input, char **output, const char *arg);
    
#ifdef __cplusplus
}
#endif

#endif
