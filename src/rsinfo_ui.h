#ifndef rstools_rsinfo_ui_h
#define rstools_rsinfo_ui_h

#include <stdio.h>
#include <strings.h>
#include <regex.h>
#include <glib.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *inputpath;
    char *dicompath;

    char *infoKey;

    char **modArgs;
    
    BOOL showComments;
    BOOL showInfo;

    rsNiftiFile *input;
    FILE *dicom;

    BOOL parametersValid;
    rsUIInterface *interface;

} rsInfoParameters;

rsInfoParameters* rsInfoParseParams(int argc, char * argv[]);
rsInfoParameters* rsInfoInitParameters();
void rsInfoBuildInterface(rsInfoParameters *p);
void rsInfoFreeParams(rsInfoParameters* p);

#ifdef __cplusplus
}
#endif

#endif
