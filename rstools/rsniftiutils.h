#include <stdio.h>
#include <strings.h>

#include <nifti1.h>
#include <fslio.h>

#if !defined(__NIFTIUTILS_H)
#define __NIFTIUTILS_H

#define BOOL    int
#define FALSE   0
#define TRUE    1

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    unsigned int x, y, z;
} Point3D;

Point3D MakePoint3D(unsigned int x, unsigned int y, unsigned int z);
double*** ResampleVolume(double ***oldVolume, int oldX, int oldY, int oldZ, int newX, int newY, int newZ);

void convertScaledDoubleToBuffer(int datatype, void *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_UINT8(  THIS_UINT8 *outbuf,   double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_INT8(   THIS_INT8 *outbuf,    double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_UINT16( THIS_UINT16 *outbuf,  double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_INT16(  THIS_INT16 *outbuf,   double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_UINT64( THIS_UINT64 *outbuf,  double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_INT64(  THIS_INT64 *outbuf,   double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_UINT32( THIS_UINT32 *outbuf,  double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_INT32(  THIS_INT32 *outbuf,   double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_FLOAT32(THIS_FLOAT32 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_FLOAT64(THIS_FLOAT64 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh);

#ifdef __cplusplus
}
#endif
    
#endif