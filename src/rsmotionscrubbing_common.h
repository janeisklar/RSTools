#ifndef rstools_rsmotionscrubbing_common_h
#define rstools_rsmotionscrubbing_common_h

#include "rsmotionscrubbing_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsMotionScrubbingInit(rsMotionScrubbingParameters *p);
void rsMotionScrubbingRun(rsMotionScrubbingParameters *p);
void rsMotionScrubbingDestroy(rsMotionScrubbingParameters *p);

void rsComputeValueRange(rsNiftiFile *file, Point3D *maskPoints, unsigned long nMaskPoints, double *min, double *max);
void rsComputeFramewiseDisplacement(double *fd, double **rp, int length);
void rsComputeDVARs(double *dvars, const double min, const double max, const rsNiftiFile *file, Point3D *maskPoints, unsigned long nMaskPoints);
void rsSaveIndexVector(char *path, BOOL *vector, const int length);
void rsSaveDoubleVector(char *path, double *vector, const int length);

inline double deg2mm(const double dist, const double deg) {
    return M_PI * dist * (deg/180.0);
}

inline int rsMax(const int a, const int b) {
  return a > b ? a : b;
}

inline int rsMin(const int a, const int b) {
  return a < b ? a : b;
}

#ifdef __cplusplus
}
#endif

#endif
