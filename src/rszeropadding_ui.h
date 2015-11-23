#ifndef rstools_zeropadding_ui_h
#define rstools_zeropadding_ui_h

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
    int padding[6];
    double paddingValue;
    BOOL mirroredPadding;

    char *callString;

    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *output;

    short newDim[3];

    rsUIInterface *interface;

} rsZeropaddingParameters;

rsZeropaddingParameters *rsZeropaddingParseParams(int argc, char * argv[]);
rsZeropaddingParameters *rsZeropaddingInitParameters();
void rsZeropaddingBuildInterface(rsZeropaddingParameters *p);
void rsZeropaddingFreeParams(rsZeropaddingParameters *p);

#ifdef __cplusplus
}
#endif

#endif
