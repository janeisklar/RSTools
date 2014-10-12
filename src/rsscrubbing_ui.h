#ifndef rstools_scrubbing_ui_h
#define rstools_scrubbing_ui_h

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
    char *flaggedpath;
    char *callString;

    BOOL *flaggedFrames;

    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *output;

    rsUIInterface *interface;

} rsScrubbingParameters;

rsScrubbingParameters *rsScrubbingParseParams(int argc, char * argv[]);
rsScrubbingParameters *rsScrubbingInitParameters();
void rsScrubbingBuildInterface(rsScrubbingParameters *p);
void rsScrubbingFreeParams(rsScrubbingParameters *p);

#ifdef __cplusplus
}
#endif

#endif
