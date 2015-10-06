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
    char *dicompath;
    char *outputpath;
    char *orientation;
    char *callString;
    char *phaseencdir;

    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *output;
    FILE *dicom;
	
    rsUIInterface *interface;

} rsOrientationParameters;

rsOrientationParameters *rsOrientationParseParams(int argc, char * argv[]);
rsOrientationParameters *rsOrientationInitParameters();
void rsOrientationBuildInterface(rsOrientationParameters *p);
void rsOrientationFreeParams(rsOrientationParameters *p);

#ifdef __cplusplus
}
#endif

#endif
