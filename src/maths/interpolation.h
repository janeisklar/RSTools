#ifndef __INTERPOLATION_MATHUTILS_H
#define __INTERPOLATION_MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

double rsInterpolationLanczosKernel(const double x, const int order);
void rsInterpolationLanczosConvolve(double* signalOut, const double* signalIn, const int nVolsIn, const int nVolsOut, const int order, const double scaling);
double rsTriLinearDistInterpolation(double***data, short xh, short yh, short zh, double voxSize[3], FloatPoint3D *point);
double rsInterpolation3DLanczosInterpolation(double***data, short xh, short yh, short zh, double voxSize[3], FloatPoint3D *point);
double rsInterpolationTriLinearInterpolation(double***data, short xh, short yh, short zh, double voxSize[3], FloatPoint3D *point);

#ifdef __cplusplus
}
#endif

#endif
