#ifndef rstools_smoothing_ui_h
#define rstools_smoothing_ui_h

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

    double kernelSizeFWHM;
    double kernelSigma;

    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *output;

    int threads;
	
    rsUIInterface *interface;
    rsReportProgressCallback *progressCallback;

} rsSmoothingParameters;

rsSmoothingParameters *rsSmoothingParseParams(int argc, char * argv[]);
rsSmoothingParameters *rsSmoothingInitParameters();
void rsSmoothingBuildInterface(rsSmoothingParameters *p);
void rsSmoothingFreeParams(rsSmoothingParameters *p);

#ifdef __cplusplus
}
#endif

#endif
