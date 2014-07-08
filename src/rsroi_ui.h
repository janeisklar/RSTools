#ifndef rstools_rsroi_ui_h
#define rstools_rsroi_ui_h

#include <stdio.h>
#include <nifti1.h>
#include <glib.h>
#include <fslio.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *inputpath;
    char *maskpath;
    char *callString;

    double sphereradius;
    FloatPoint3D *center;
    FloatPoint3D *cubeDim;
    BOOL keepVolume;
    long nSamples;
    double roiValue;
    BOOL useImageSpaceCoordinates;

    BOOL verbose;
    BOOL parametersValid;

    GOptionContext *context;

    rsNiftiFile *input;
    rsNiftiFile *mask;

    BOOL inputEqualsOutput;

} rsRoiParameters;

rsRoiParameters *rsRoiParseParams(int argc, char * argv[]);
rsRoiParameters *rsRoiInitParameters();
void rsRoiFreeParams(rsRoiParameters *p);
void rsRoiPrintHelp(rsRoiParameters *p);

gboolean rsRoiParsePoint(const gchar *option_name, const gchar *value, gpointer data, GError **error);

#ifdef __cplusplus
}
#endif

#endif
