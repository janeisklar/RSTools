#ifndef rstools_motionscrubbing_ui_h
#define rstools_motionscrubbing_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RSTOOLS_REALIGNMENT_PARAMETER_FORMAT_SPM 1
#define RSTOOLS_REALIGNMENT_PARAMETER_FORMAT_FSL 2

typedef struct {
    char *inputpath;
    char *maskpath;
    char *outputpath;
    char *realignmentpath;
    char *fdpath;
    char *dvarspath;
    char *flaggedpath;
    char *callString;

    double fdthreshold;
    double dvarsthreshold;

    long rpColumns;
    long rpEntries;
    double **rp;
    Point3D *maskPoints;
    unsigned long nMaskPoints;

    short rpformat;

    BOOL useModal;
    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *output;
    double ***mask;

    rsUIInterface *interface;

    int threads;

} rsMotionScrubbingParameters;

rsMotionScrubbingParameters *rsMotionScrubbingParseParams(int argc, char * argv[]);
rsMotionScrubbingParameters *rsMotionScrubbingInitParameters();
void rsMotionScrubbingBuildInterface(rsMotionScrubbingParameters *p);
void rsMotionScrubbingFreeParams(rsMotionScrubbingParameters *p);

gboolean rsMotionScrubbingParseRPFormat(const gchar *option_name, const gchar *value, gpointer data, GError **error);

#ifdef __cplusplus
}
#endif

#endif
