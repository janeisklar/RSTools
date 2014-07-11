#include <stdio.h>
#include <strings.h>
#include <ctype.h>
#include <time.h>

#include <omp.h>

#include <nifti1.h>
#include <fslio.h>
#include <gsl/gsl_statistics_double.h>
#include "src/rscommon.h"
#include "src/utils/rsmemory.h"

#if !defined(__NIFTIUTILS_H)
#define __NIFTIUTILS_H

#ifdef __cplusplus
extern "C" {
#endif
    
#if !defined(BOOL)
#define BOOL    int
#endif

#if !defined(FALSE)
#define FALSE   0
#endif

#if !defined(TRUE)
#define TRUE    1
#endif
    
#define RSIOERR(x) { fprintf(stderr,"Error:: %s\n",(x)); fflush(stderr); exit(EXIT_FAILURE); }

#define RSNIFTI_OPEN_NONE     1   // file will be opened only and the header info read
#define RSNIFTI_OPEN_ALLOC    2   // additionally the space required to load/save the file is allocated
#define RSNIFTI_OPEN_READ     3   // additionally the file is also read into the allocated space
#define RSNIFTI_CLONE_POINTER 4   // don't allocate any space, just re-use the pointer from the original volume for the clone
#define RSNIFTI_CLONE_MEMORY  5   // allocate the required space and copy it the image content from the original volume over to the clone

#define RSNIFTI_CLONE_AS_INPUT 0  // use the same setting as in the volume that is to be cloned

typedef struct {
    unsigned int x, y, z;
} Point3D;
    
typedef struct {
    float x, y, z;
} FloatPoint3D;

typedef struct {
    char*   path;
    BOOL    readable;
    FSLIO*  fslio;
    short   xDim;
    short   yDim;
    short   zDim;
    short   vDim;
    short   dt;
    size_t  pixtype;
    float   inter;
    float   slope;
    void*   data;

    short   intent_code;
    float   intent_p1;
    float   intent_p2;
    float   intent_p3;
} rsNiftiFile;

typedef struct {
    char*         originalPath;
    char*         resampledPath;
    unsigned long nPoints;
    Point3D*      maskPoints;
    double***     resampledMask;
    BOOL          readable;
} rsMask;

Point3D*      rsMakePoint3D(unsigned int x, unsigned int y, unsigned int z);
FloatPoint3D* rsMakeFloatPoint3D(float x, float y, float z);
BOOL          rsPointInVolume(const Point3D *p, const int xh, const int yh, const int zh);
Point3D*      rsReadMask(char *path, unsigned short newX, unsigned short newY, unsigned short newZ, unsigned long *nPoints, char *resampledMaskPath, FSLIO *maskPrototype, double ***resampledMaskReturn);
rsMask*       rsMaskInit(char *path);
void          rsMaskLoad(rsMask *mask, rsNiftiFile *resamplingPrototype);
void          rsMaskFree(rsMask *mask);
long          rsCleanMaskFromNaNs(rsMask *mask, rsNiftiFile *input);
double***     rsResampleVolume(double ***oldVolume, int oldX, int oldY, int oldZ, int newX, int newY, int newZ);
size_t        rsWriteTimeSeries(FSLIO *fslio, const void *buffer, short xVox, short yVox, short zVox, int nvols);
BOOL          rsExtractTimecourseFromBuffer(const short datatype, double *timecourse, const void *buffer, const float slope, const float inter, const Point3D *p, const int xh, const int yh, const int zh, const int th);
BOOL          rsExtractPointsFromBuffer(const short datatype, double *data, const void *buffer, const float slope, const float inter, const Point3D *points, const unsigned long nPoints, const int t, const int xh, const int yh, const int zh, const int th);
BOOL          rsWriteTimecourseToBuffer(const short datatype, const double *timecourse, void *buffer, const float slope, const float inter, const Point3D *p, const int xh, const int yh, const int zh, const int th);
BOOL          rsResetBufferToValue(const short datatype, void *buffer, const float slope, const float inter, const int xh, const int yh, const int zh, const int th, const double value);
BOOL          rsWriteVolumeToBuffer(const short datatype, double *data, const void *buffer, const float slope, const float inter, const int t, const int xh, const int yh, const int zh);
BOOL          rsExtractVolumeFromBuffer(const short datatype, double *data, const void *buffer, const float slope, const float inter, const int t, const int xh, const int yh, const int zh);
size_t        rsWordLength(const short datatype);
size_t        rsVolumeOffset(const int xh, const int yh, const int zh, const int t);
size_t        rsVolumeLength(const int xh, const int yh, const int zh);
size_t        rsGetBufferSize(const int xh, const int yh, const int zh, const int th, const int nifti_datatype);

// the offset of a voxel within a volume(in words)
inline size_t rsVoxelOffset(const size_t x, const size_t y, const size_t z, const size_t xh, const size_t yh) {
    return z * (xh * yh) + (y * xh) + x;
}

// the overall offset of a voxel at a certain time frame within a 4D-dataset
inline size_t rsOverallVoxelOffset(const size_t x, const size_t y, const size_t z, const size_t t, const size_t xh, const size_t yh, const size_t zh) {
    return t * (xh * yh * zh) + z * (xh * yh) + (y * xh) + x;
}
    
BOOL convertScaledDoubleToBuffer(const short datatype, void *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_UINT8(  THIS_UINT8 *outbuf,   double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_INT8(   THIS_INT8 *outbuf,    double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_UINT16( THIS_UINT16 *outbuf,  double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_INT16(  THIS_INT16 *outbuf,   double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_UINT64( THIS_UINT64 *outbuf,  double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_INT64(  THIS_INT64 *outbuf,   double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_UINT32( THIS_UINT32 *outbuf,  double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_INT32(  THIS_INT32 *outbuf,   double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_FLOAT32(THIS_FLOAT32 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh);
void convertScaledDoubleToBuffer_FLOAT64(THIS_FLOAT64 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh);

char *rsLeftTrimString(char *s);
char *rsRightTrimString(char *s);
char *rsTrimString(char *s);
char *rsMergeStringArray(int argc, char * argv[]);
char *rsReadCommentFile(char *path);

void rsSetThreadsNum(const unsigned int threads);
unsigned int rsGetThreadsNum();
    
rsNiftiFile *rsOpenNiftiFile(const char* path, const unsigned int mode);
rsNiftiFile *rsInitNiftiFile(void);
void         rsFreeNiftiFile(rsNiftiFile* f);
void         rsCloseNiftiFileAndFree(rsNiftiFile* f);
rsNiftiFile *rsCloneNiftiFile(const char* path, const rsNiftiFile* f, const unsigned int mode, const int vDim);
void        rsCloseNiftiFile(rsNiftiFile* f, BOOL keepData);
void        rsWriteNiftiHeader(FSLIO *fslio, char* description);

#define rsWriteTimecourseToRSNiftiFileBuffer(nifti,convdata,point);                 rsWriteTimecourseToBuffer(nifti->dt, convdata, nifti->data, nifti->slope, nifti->inter, point, nifti->xDim, nifti->yDim, nifti->zDim, nifti->vDim);
#define rsExtractTimecourseFromRSNiftiFileBuffer(nifti,convdata,point);             rsExtractTimecourseFromBuffer(nifti->dt, convdata, nifti->data, nifti->slope, nifti->inter, point, nifti->xDim, nifti->yDim, nifti->zDim, nifti->vDim);
#define rsExtractVolumeFromRSNiftiFileBuffer(nifti,convdata,t);                     rsExtractVolumeFromBuffer(nifti->dt, convdata, nifti->data, nifti->slope, nifti->inter, t, nifti->xDim, nifti->yDim, nifti->zDim);
#define rsWriteVolumeToRSNiftiFileBuffer(nifti,convdata,t);                         rsWriteVolumeToBuffer(nifti->dt, convdata, nifti->data, nifti->slope, nifti->inter, t, nifti->xDim, nifti->yDim, nifti->zDim);
#define rsResetRSNiftiFileBufferToValue(nifti,value);                               rsResetBufferToValue(nifti->dt, nifti->data, nifti->slope, nifti->inter, nifti->xDim, nifti->yDim, nifti->zDim, nifti->vDim, value);
#define rsExtractPointsFromRSNiftiFileBuffer(Tnifti,Tconvdata,Tpoints,TnPoints,Tt); rsExtractPointsFromBuffer(Tnifti->dt, Tconvdata, Tnifti->data, Tnifti->slope, Tnifti->inter, Tpoints, TnPoints, Tt, Tnifti->xDim, Tnifti->yDim, Tnifti->zDim, Tnifti->vDim);


#ifdef __cplusplus
}
#endif
    
#endif
