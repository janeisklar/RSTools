#ifndef __INTERPOLATION_MATHUTILS_H
#define __INTERPOLATION_MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

double rsInterpolationLanczosKernel(const double x, const int order);
void rsInterpolationLanczosConvolve(double* signalOut, const double* signalIn, const int nVolsIn, const int nVolsOut, const int order, const double scaling);
double rsInterpolationLinear(double x, double x1, double x2, double q00, double q01);
double rsInterpolationBiLinear(double x, double y, double q11, double q12, double q21, double q22, double x1, double x2, double y1, double y2);
double rsInterpolationTriLinear(double x, double y, double z, double q000, double q001, double q010, double q011, double q100, double q101, double q110, double q111, double x1, double x2, double y1, double y2, double z1, double z2);
double rsTriLinearDistInterpolation(double***data, short xh, short yh, short zh, FloatPoint3D *point);

#ifdef __cplusplus
}
#endif

#endif
