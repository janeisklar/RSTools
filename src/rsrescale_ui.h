#ifndef rstools_rescale_ui_h
#define rstools_rescale_ui_h

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

    double scale[3];
    double appliedScale[3];
    float oldSpacing[3];
    float newSpacing[3];
    short newDim[3];
    float TR;

    mat44 outputWorldMatrix;
    mat44 invInputWorldMatrix;

    BOOL linearInterpolation;
    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *output;

    int threads;
	
    rsUIInterface *interface;
    rsReportProgressCallback *progressCallback;

} rsRescaleParameters;

rsRescaleParameters *rsRescaleParseParams(int argc, char * argv[]);
rsRescaleParameters *rsRescaleInitParameters();
void rsRescaleBuildInterface(rsRescaleParameters *p);
void rsRescaleFreeParams(rsRescaleParameters *p);

#ifdef __cplusplus
}
#endif

#endif
