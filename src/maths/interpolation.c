#include <nifti/rsniftiutils.h>
#include <src/nifti/rsniftiutils.h>
#include "interpolation.h"
#include "math.h"


double rsInterpolationLanczosKernel(const double x, const int order)
{
    if (x == 0) {
        return 1.0;
    } else if (order > fabs(x) && fabs(x) > 0.0) {
	    return order * sin(M_PI * x) * sin(M_PI * x / order) / (pow(M_PI, 2.0) * pow(x, 2.0));
    } else {
        return 0.0;
    }
}

void rsInterpolationLanczosConvolve(double* signalOut, const double* signalIn, const int nVolsIn, const int nVolsOut, const int order, const double scaling)
{
    int j = 0;
    for (double x=order-1; x<=nVolsIn-order+1; x+=1/scaling) {

        const int xi = (int)floor(x);
        double S=0;
        
        for (int i=xi-order+1; i<xi+order; i++) {
            S += signalIn[i]*rsInterpolationLanczosKernel(xi-i, order);
        }
        signalOut[j++] = S;
    }
}

double rsInterpolationLinear(double x, double x1, double x2, double q00, double q01)
{
    return ((x2 - x) / (x2 - x1)) * q00 + ((x - x1) / (x2 - x1)) * q01;
}

double rsInterpolationBiLinear(double x, double y, double q11, double q12, double q21, double q22, double x1, double x2, double y1, double y2) {
    const double r1 = rsInterpolationLinear(x, x1, x2, q11, q21);
    const double r2 = rsInterpolationLinear(x, x1, x2, q12, q22);

    return rsInterpolationLinear(y, y1, y2, r1, r2);
}

double rsInterpolationTriLinear(double x, double y, double z, double q000, double q001, double q010, double q011, double q100, double q101, double q110, double q111, double x1, double x2, double y1, double y2, double z1, double z2)
{
    const double x00 = rsInterpolationLinear(x, x1, x2, q000, q100);
    const double x10 = rsInterpolationLinear(x, x1, x2, q010, q110);
    const double x01 = rsInterpolationLinear(x, x1, x2, q001, q101);
    const double x11 = rsInterpolationLinear(x, x1, x2, q011, q111);
    const double r0 = rsInterpolationLinear(y, y1, y2, x00, x01);
    const double r1 = rsInterpolationLinear(y, y1, y2, x10, x11);

    return rsInterpolationLinear(z, z1, z2, r0, r1);
}

double rsInterpolationTriLinearInterpolation(double***data, short xh, short yh, short zh, double voxSize[3], FloatPoint3D *point) {
    SignedPoint3D *a = rsMakeSignedPoint3D(floor(point->x), floor(point->y), floor(point->z));
    SignedPoint3D *b = rsMakeSignedPoint3D(a->x+1, a->y+1, a->z+1);

    // construct 8 points surrounding the point to be interpolated
    SignedPoint3D* points[8] = {
        rsMakeSignedPoint3D(a->x, a->y, a->z), // q000
        rsMakeSignedPoint3D(b->x, a->y, a->z), // q001
        rsMakeSignedPoint3D(a->x, b->y, a->z), // q010
        rsMakeSignedPoint3D(b->x, b->y, a->z), // q011
        rsMakeSignedPoint3D(a->x, a->y, b->z), // q100
        rsMakeSignedPoint3D(b->x, a->y, b->z), // q101
        rsMakeSignedPoint3D(a->x, b->y, b->z), // q110
        rsMakeSignedPoint3D(b->x, b->y, b->z)  // q111
    };

    // return NaN if point exceeds the volume
    for (short i=0; i<8; i++) {
        if (!rsSignedPointInVolume(points[i], xh, yh, zh)) {
            return log(-1.0);
        }
    }

    // compute interpolation between the specified points
    return rsInterpolationTriLinear(
        point->x,
        point->y,
        point->z,
        data[points[0]->z][points[0]->y][points[0]->x], // q000
        data[points[1]->z][points[1]->y][points[1]->x], // q001
        data[points[2]->z][points[2]->y][points[2]->x], // q010
        data[points[3]->z][points[3]->y][points[3]->x], // q011
        data[points[4]->z][points[4]->y][points[4]->x], // q100
        data[points[5]->z][points[5]->y][points[5]->x], // q101
        data[points[6]->z][points[6]->y][points[6]->x], // q110
        data[points[7]->z][points[7]->y][points[7]->x], // q111
        a->x,
        b->x,
        a->y,
        b->y,
        a->z,
        b->z
    );
}

double rsInterpolation3DLanczosInterpolation(double***data, short xh, short yh, short zh, double voxSize[3], FloatPoint3D *point) {
    const int order = 3;

    SignedPoint3D *convStart = rsMakeSignedPoint3D(
        floor(point->x)-order+1,
        floor(point->y)-order+1,
        floor(point->z)-order+1
    );
    SignedPoint3D *convEnd = rsMakeSignedPoint3D(
        floor(point->x)+order,
        floor(point->y)+order,
        floor(point->z)+order
    );

    double result = 0.0;

    // return a linear interpolated value if we exceed the bounds of this interpolation method
    if (!rsSignedPointInVolume(convStart, xh, yh, zh) || !rsSignedPointInVolume(convEnd, xh, yh, zh)) {
        rsFree(convStart);
        rsFree(convEnd);
        return rsTriLinearDistInterpolation(data, xh, yh, zh, voxSize, point);
    }

    for (short i=convStart->x; i <= convEnd->x; i++) {
        for (short j=convStart->y; j <= convEnd->y; j++) {
            for (short k=convStart->z; k <= convEnd->z; k++) {
                const double S = data[k][j][i];

                // use linear interpolation when dealing with NaN/Inf values
                if(!isfinite(S)) {
                    return rsTriLinearDistInterpolation(data, xh, yh, zh, voxSize, point);
                }

                const double Li = rsInterpolationLanczosKernel(point->x-i, order);
                const double Lj = rsInterpolationLanczosKernel(point->y-j, order);
                const double Lk = rsInterpolationLanczosKernel(point->z-k, order);

                result += S * Li * Lj * Lk;
            }
        }
    }

    rsFree(convStart);
    rsFree(convEnd);

    return result;
}

double rsTriLinearDistInterpolation(double***data, short xh, short yh, short zh, double voxSize[3], FloatPoint3D *point) {
    SignedPoint3D *a = rsMakeSignedPoint3D(floor(point->x), floor(point->y), floor(point->z));
    SignedPoint3D *b = rsMakeSignedPoint3D(a->x+1, a->y+1, a->z+1);

    // construct 8 points surrounding the point to be interpolated
    SignedPoint3D* points[8] = {
        rsMakeSignedPoint3D(a->x, a->y, a->z), // q000
        rsMakeSignedPoint3D(b->x, a->y, a->z), // q001
        rsMakeSignedPoint3D(a->x, b->y, a->z), // q010
        rsMakeSignedPoint3D(b->x, b->y, a->z), // q011
        rsMakeSignedPoint3D(a->x, a->y, b->z), // q100
        rsMakeSignedPoint3D(b->x, a->y, b->z), // q101
        rsMakeSignedPoint3D(a->x, b->y, b->z), // q110
        rsMakeSignedPoint3D(b->x, b->y, b->z)  // q111
    };

    // determine which points lie within the volume and compute their distances
    double summedInvDistance=0.0;
    double weightedInvSum = 0.0;
    short nValidPoints = 0;
    for (short i=0; i<8; i++) {
        BOOL validPoint = rsSignedPointInVolume(points[i], xh, yh, zh);
        if (validPoint) {

            const double S = data[points[i]->z][points[i]->y][points[i]->x];

            // ignore NaN and Inf values
            if (!isfinite(S)) {
                continue;
            }

            // inv dist in each direction
            float invDist[3] = {
                1.0 - fabs(point->x - points[i]->x),
                1.0 - fabs(point->y - points[i]->y),
                1.0 - fabs(point->z - points[i]->z)
            };

            // scale by the corresponding voxel size
            for (short d=0; d<3; d++) {
                invDist[d] = invDist[d] / voxSize[d];
            }

            // euclidean dist sqrt((1-dx)^2 + (1-dy)^2 + (1-dz)^2)
            const double invDistance = pow(
                (
                    pow(invDist[0], 2.0) +
                    pow(invDist[1], 2.0) +
                    pow(invDist[2], 2.0)
                ),
                0.5
            );

            summedInvDistance += invDistance;
            weightedInvSum += invDistance * S;

            nValidPoints++;
        }
    }

    // clear data
    for (short i=0; i<8; i++) {
        rsFree(points[i]);
    }
    rsFree(a);
    rsFree(b);

    // terminate with NaN if no points lie withing the volume
    if (nValidPoints < 1) {
        return log(-1.0);
    }

    // normalize weighted sum and return it
    return weightedInvSum / summedInvDistance;
}
