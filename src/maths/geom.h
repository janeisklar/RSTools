#include <stdio.h>
#include <math.h>

#include "src/nifti/rsniftiutils.h"

#if !defined(__GEOM_MATHUTILS_H)
#define __GEOM_MATHUTILS_H

#ifdef __cplusplus
extern "C" {
#endif 

BOOL rsVoxelInSphere(FloatPoint3D point, FloatPoint3D center, double radius);
BOOL rsVoxelInCube(FloatPoint3D point, FloatPoint3D center, FloatPoint3D dim);
double rsDistance(FloatPoint3D A, FloatPoint3D B);

#endif
