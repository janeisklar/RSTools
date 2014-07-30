#ifndef rstools_rsroi_ui_h
#define rstools_rsroi_ui_h

#include <stdio.h>
#include <fslio.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"
#include "src/utils/rsui.h"

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
    
    rsUIInterface *interface;

    rsNiftiFile *input;
    rsNiftiFile *mask;

    BOOL inputEqualsOutput;

} rsRoiParameters;

rsRoiParameters *rsRoiParseParams(int argc, char * argv[]);
rsRoiParameters *rsRoiInitParameters();
void rsRoiBuildInterface(rsRoiParameters *p);
void rsRoiFreeParams(rsRoiParameters *p);

gboolean rsRoiParsePoint(const gchar *option_name, const gchar *value, gpointer data, GError **error);

#ifdef __cplusplus
}
#endif

#endif
