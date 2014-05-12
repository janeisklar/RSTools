#include "geom.h"

double rsDistance(FloatPoint3D A, FloatPoint3D B)
{
    return pow(
        (
            pow(A.x-B.x,2.0) +
            pow(A.y-B.y,2.0) +
            pow(A.z-B.z,2.0)
        ),
        1.0/2.0
    );
}

BOOL rsVoxelInSphere(FloatPoint3D point, FloatPoint3D center, double radius)
{
    return rsDistance(point, center) <= radius;
}

BOOL rsVoxelInCube(FloatPoint3D point, FloatPoint3D center, FloatPoint3D dim)
{
    return fabs(point.x-center.x) <= (dim.x / 2.0) &&
           fabs(point.y-center.y) <= (dim.y / 2.0) &&
           fabs(point.z-center.z) <= (dim.z / 2.0);
}