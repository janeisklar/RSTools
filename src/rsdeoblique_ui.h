#ifndef rstools_deoblique_ui_h
#define rstools_deoblique_ui_h

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
    char *transformationpath;

    char *callString;

    BOOL verbose;
    BOOL parametersValid;
    BOOL retainVoxelSizes;

    rsNiftiFile *input;
    rsNiftiFile *output;
    FILE *transform;

    short dimsOut[3];
    double pixDimIn[3];
    double pixDimOut[3];
    float TR;
    short padding;

    rsUIInterface *interface;

} rsDeobliqueParameters;

rsDeobliqueParameters *rsDeobliqueParseParams(int argc, char * argv[]);
rsDeobliqueParameters *rsDeobliqueInitParameters();
void rsDeobliqueBuildInterface(rsDeobliqueParameters *p);
void rsDeobliqueFreeParams(rsDeobliqueParameters *p);

#ifdef __cplusplus
}
#endif

#endif
