#ifndef rstools_rstimecourse_ui_h
#define rstools_rstimecourse_ui_h

#include <stdio.h>
#include <strings.h>
#include <dlfcn.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#define RSTOOLS_TIMECOURSE_ALGORITHM_MEAN   1
#define RSTOOLS_TIMECOURSE_ALGORITHM_SPCA   2
#define RSTOOLS_TIMECOURSE_ALGORITHM_TPCA   3
#define RSTOOLS_TIMECOURSE_ALGORITHM_CSP    4
#define RSTOOLS_TIMECOURSE_ALGORITHM_STDDEV 5

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    
    char *inputpath;
    char *outputpath;
    char *maskpath;
    char *mask2path;
    char *savemaskpath;
    char *eigenvaluespath;
    char *spatialmappath;
    
    short algorithm;
    
    BOOL useStandardScores;
    
    float minVariance;
    int nComponents;
    
    rsNiftiFile *input;
    FILE *output;
    rsMask *mask;
    rsMask *mask2;
    Point3D *point;
    
    BOOL parametersValid;
    BOOL verbose;
    rsUIInterface *interface;
    char *callString;
    int threads;

} rsTimecourseParameters;

rsTimecourseParameters* rsTimecourseParseParams(int argc, char * argv[]);
rsTimecourseParameters* rsTimecourseInitParameters();
void rsTimecourseBuildInterface(rsTimecourseParameters *p);
void rsTimecourseFreeParams(rsTimecourseParameters* p);

gboolean rsTimecourseParseAlgorithm(const gchar *option_name, const gchar *value, gpointer data, GError **error);
gboolean rsTimecourseParsePoint(const gchar *option_name, const gchar *value, gpointer data, GError **error);

#ifdef __cplusplus
}
#endif

#endif
