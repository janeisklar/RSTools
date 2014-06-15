#ifndef rstools_rstimecourse_common_h
#define rstools_rstimecourse_common_h

#include "rstimecourse_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsTimecourseInit(rsTimecourseParameters* p);
void rsTimecourseRun(rsTimecourseParameters *p);
void rsTimecourseDestroy(rsTimecourseParameters* p);

void rsTimecourseRunSingleVoxelExtraction(rsTimecourseParameters *p);
void rsTimecourseRunMeanOrStdDev(rsTimecourseParameters *p);
void rsTimecourseRunPCA(rsTimecourseParameters *p);
void rsTimecourseRunCSP(rsTimecourseParameters *p);

void rsWriteSpatialMap(char *file, const rsNiftiFile *reference, Point3D *points, gsl_matrix *maps);

#ifdef __cplusplus
}
#endif

#endif
