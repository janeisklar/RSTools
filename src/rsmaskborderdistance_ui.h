#ifndef rstools_maskborderdistance_ui_h
#define  rstools_maskborderdistance_ui_h

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
    
    BOOL verbose;
    BOOL parametersValid;
    
    rsNiftiFile *input;
    rsNiftiFile *output;

    rsUIInterface *interface;

    int threads;

    rsReportProgressCallback *progressCallback;
    
} rsMaskBorderDistanceParameters;

rsMaskBorderDistanceParameters *rsMaskBorderDistanceParseParams(int argc, char * argv[]);
rsMaskBorderDistanceParameters *rsMaskBorderDistanceInitParameters();
void rsMaskBorderDistanceFreeParams(rsMaskBorderDistanceParameters *p);
void rsMaskBorderDistanceBuildInterface(rsMaskBorderDistanceParameters *p);
    
#ifdef __cplusplus
}
#endif

#endif
