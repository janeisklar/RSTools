#include "rsniftiutils.h"

static unsigned int rsThreadsNum = 1;

void rsSetThreadsNum(const unsigned int threads)
{
	rsThreadsNum = threads;
}

unsigned int rsGetThreadsNum()
{
	return rsThreadsNum;
}

/*
 * Creates a new 3D point
 */
Point3D MakePoint3D(unsigned int x, unsigned int y, unsigned int z)
{
    Point3D a;
    a.x = x;
    a.y = y;
    a.z = z;
    return a;
}

FloatPoint3D MakeFloatPoint3D(float x, float y, float z)
{
    FloatPoint3D a;
    a.x = x;
    a.y = y;
    a.z = z;
    return a;
}

/*
 * Reads in a binary mask in terms of a nifti file and
 * rescales it to the specified dimensions. It then
 * returns an array with all the points in the mask
 * that have a value higher than 0.
 * The resulting mask after rescaling can be saved as
 * a nifti by supplying an additional path. The headers
 * for this file will be cloned from the given 
 * prototype. This can be used to ensure that the used
 * mask is appropriate for the orientation of the file
 * that it is later applied to.
 */
Point3D* ReadMask(char *path, unsigned short newX, unsigned short newY, unsigned short newZ, unsigned long *nPoints, char *resampledMaskPath, FSLIO *maskPrototype, double ***resampledMaskReturn)
{
    FSLIO *fslio;
	void *buffer;
	unsigned long buffsize;
    
    short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    
    int pointCount = 0;
    
    /* Open mask */
    fslio = FslOpen(path, "rb");
    if (fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",path);
        return NULL;
    }
    
    /* Read out dimensions */
    FslGetDim(fslio, &xDim, &yDim, &zDim, &vDim);
    
    if (fslio->niftiptr->scl_slope != 0) {
        slope = fslio->niftiptr->scl_slope;
        inter = fslio->niftiptr->scl_inter;
    }
    
    /* Determine datatype */
	dt = FslGetDataType(fslio, &pixtype);
    
    /* Init buffer */
    buffsize = (unsigned long)xDim * (unsigned long)yDim * (unsigned long)zDim * (unsigned long)(dt/8);
    buffer   = malloc(buffsize);
    
    /* Read in first volume */
    if (!FslReadVolumes(fslio, buffer, 1)) {
        free(buffer);
        fprintf(stderr, "\nError - reading data in %s\n", path);
        FslClose(fslio);
        return NULL;
    }
    
    double ***mask = FslGetVolumeAsScaledDouble(fslio, 0);
    
    /* Resample mask to have the same scaling as the input volume */
    double ***resampledMask = ResampleVolume(mask, xDim, yDim, zDim, newX, newY, newZ);
    
    free(mask[0][0]);
    free(mask[0]);
    free(mask);
    
    /* Count how many points we'll get */
    *nPoints = (unsigned long)0L;
	unsigned int x=0,y=0,z=0;
    for (x=0; x<newX; x=x+1) {
        for (y=0; y<newY; y=y+1) {
            for (z=0; z<newZ; z=z+1) {
                if ( resampledMask[z][y][x] > 0.01 ) {
                    *nPoints = *nPoints + ((unsigned long)1L);
                }
                
                if ( resampledMaskReturn != NULL ) {
                    resampledMaskReturn[z][y][x] = resampledMask[z][y][x];
                }
            }
        }
    }
    
    /* Initialize result array */
    Point3D* points = malloc(*nPoints*sizeof(Point3D));
    
    /* Create array with all points that are in the mask */
    int i=0;
    for (x=0; x<newX; x=x+1) {
        for (y=0; y<newY; y=y+1) {
            for (z=0; z<newZ; z=z+1) {
                if ( resampledMask[z][y][x] > 0.01 ) {
                    points[i] = MakePoint3D(x,y,z);
                    i = i+1;
                }
            }
        }
    }
    
    /* Save mask */
    FSLIO *fslioResampled = NULL;
    
    if ( resampledMaskPath != NULL ) {
        fslioResampled = FslOpen(resampledMaskPath, "wb");
    }
    if (fslioResampled == NULL) {
        if ( resampledMaskPath != NULL ) {
            fprintf(stderr, "\nWarning, could not open %s for writing.\n",resampledMaskPath);
        }
    } else {
        
        dt = FslGetDataType(maskPrototype, &pixtype);
        
        FslCloneHeader(fslioResampled, maskPrototype);
        FslSetDim(fslioResampled, newX,newY,newZ,1);
        FslSetDimensionality(fslioResampled, 3);
        FslSetDataType(fslioResampled, pixtype);
        FslWriteHeader(fslioResampled);
        
        void *maskBuffer = malloc((unsigned long)newX * (unsigned long)newY * (unsigned long)newZ * (unsigned long)(dt/8));
        
        convertScaledDoubleToBuffer(
            maskPrototype->niftiptr->datatype,
            maskBuffer,
            resampledMask[0][0],
            maskPrototype->niftiptr->scl_slope,
            maskPrototype->niftiptr->scl_inter,
            newX,
            newY,
            newZ,
            TRUE
        );
        
        FslWriteVolumes(fslioResampled, maskBuffer, 1);
        FslClose(fslioResampled);
    }
    
    FslClose(fslio);
    free(buffer);
    free(resampledMask[0][0]);
    free(resampledMask[0]);
    free(resampledMask);
    free(fslio);
    free(fslioResampled);
    
    return points;
}

/*
 * Resamples a 3D volume from a given resolution into a new resolution. 
 * The 3D volume thereby needs to be supplied as a 3D matrix of scaled
 * doubles as created by d3matrix.
 */
double*** ResampleVolume(double ***oldVolume, int oldX, int oldY, int oldZ, int newX, int newY, int newZ)
{
    double ***resampledVolume = d3matrix(newZ, newY, newX);
    
    for (int x=0; x<newX; x=x+1) {
        
        int resampledX = ((float)oldX / (float)newX) * (float)x;
        for (int y=0; y<newY; y=y+1) {
            
            int resampledY = ((float)oldY / (float)newY) * (float)y;
            for (int z=0; z<newZ; z=z+1) {
                
                int resampledZ = ((float)oldZ / (float)newZ) * (float)z;
                resampledVolume[z][y][x] = oldVolume[resampledZ][resampledY][resampledX];
            }
        }
    }
    
    return resampledVolume;
}

size_t rsWriteTimeSeries(FSLIO *fslio, const void *buffer, short xVox, short yVox, short zVox, int nvols)
{
    size_t volbytes, offset, orig_offset;
    size_t n;
    short xdim,ydim,zdim,v,wordsize;
    
    if (fslio==NULL)  fprintf(stderr, "rsWriteTimeSeries: Null pointer passed for FSLIO");
    if (fslio->niftiptr!=NULL) {
        
        FslGetDim(fslio,&xdim,&ydim,&zdim,&v);
        
        if ((xVox<0) || (xVox >=xdim)){
            fprintf(stderr, "rsWriteTimeSeries: voxel coordinate(%hd) outside valid x-range(%hd..%hd)", xVox, 0, xdim-1);
            return 0;
        }
        if ((yVox<0) || (yVox >=ydim)) {
            fprintf(stderr, "rsWriteTimeSeries: voxel coordinate(%hd) outside valid y-range(%hd..%hd)", yVox, 0, ydim-1);
            return 0;
        }
        if ((zVox<0) || (zVox >=zdim)) {
            fprintf(stderr, "rsWriteTimeSeries: voxel coordinate(%hd) outside valid z-range(%hd..%hd)", zVox, 0, zdim-1);
            return 0;
        }
        
        wordsize = fslio->niftiptr->nbyper;
        volbytes = xdim * ydim * zdim * wordsize;
        
        orig_offset = znztell(fslio->fileptr);
        
        /* go back to the beginning of the data(after the header) */
        znzseek(fslio->fileptr, fslio->niftiptr->iname_offset, SEEK_SET);
        offset = ((ydim * zVox + yVox) * xdim + xVox) * wordsize;
        znzseek(fslio->fileptr,offset,SEEK_CUR);
        
        for (n=0; n<nvols; n++) {
            if (n>0) znzseek(fslio->fileptr, volbytes - wordsize, SEEK_CUR);
/*            fprintf(stdout, "Wrote %03zd: %hd.\n", n, (short)*((char*)(buffer)+(n*wordsize))); */
            if (znzwrite((char *)buffer+(n*wordsize), 1, wordsize,fslio->fileptr) != wordsize) {
                fprintf(stderr, "rsWriteTimeSeries: failed to write values");
                return 0;
            }
            /*if (fslio->niftiptr->byteorder != nifti_short_order())
                nifti_swap_Nbytes(1,fslio->niftiptr->swapsize,
                                  (char *)buffer+(n*wordsize));*/
        }
        
        /* restore file pointer to original position */
        znzseek(fslio->fileptr,orig_offset,SEEK_SET);
        return n;
    }
    if (fslio->mincptr!=NULL) {
        fprintf(stderr,"Warning:: Minc is not yet supported\n");
    }
    return 0;
}

void rsScaleValues(double *timecourse, const int th, const double slope, const double inter) {
    for (int t=0; t<th; t=t+1) {
        timecourse[t] = timecourse[t] * slope + inter;
    }
}

/*
 * Takes in a buffer obtained from fsl and reads
 * out a timecourse for a single voxel that is
 * specified by the parameter p.
 * The variables xh, yh, zh, and th correspond
 * to the dimensions of the nifti file.
 *
 */
BOOL rsExtractTimecourseFromBuffer(const FSLIO *fslio, double *timecourse, const void *buffer, const float slope, const float inter, const Point3D p, const int xh, const int yh, const int zh, const int th) {
    
    const long voxeLOffset = p.z * (xh * yh) + (p.y * xh) + p.x;
    const long timeOffset  = xh * yh * zh;
    
    for (int t=0; t<th; t=t+1) {
        
        const long address = t * timeOffset + voxeLOffset;
        
        switch (fslio->niftiptr->datatype) {
            case NIFTI_TYPE_UINT8:
                timecourse[t] = (double) *((THIS_UINT8 *)(buffer)+address);
                break;
            case NIFTI_TYPE_INT8:
                timecourse[t] = (double) *((THIS_INT8 *)(buffer)+address);
                break;
            case NIFTI_TYPE_UINT16:
                timecourse[t] = (double) *((THIS_UINT16 *)(buffer)+address);
                break;
            case NIFTI_TYPE_INT16:
                timecourse[t] = (double) *((THIS_INT16 *)(buffer)+address);
                break;
            case NIFTI_TYPE_UINT64:
                timecourse[t] = (double) *((THIS_UINT64 *)(buffer)+address);
                break;
            case NIFTI_TYPE_INT64:
                timecourse[t] = (double) *((THIS_INT64 *)(buffer)+address);
                break;
            case NIFTI_TYPE_UINT32:
                timecourse[t] = (double) *((THIS_UINT32 *)(buffer)+address);
                break;
            case NIFTI_TYPE_INT32:
                timecourse[t] = (double) *((THIS_INT32 *)(buffer)+address);
                break;
            case NIFTI_TYPE_FLOAT32:
                timecourse[t] = (double) *((THIS_FLOAT32 *)(buffer)+address);
                break;
            case NIFTI_TYPE_FLOAT64:
                timecourse[t] = (double) *((THIS_FLOAT64 *)(buffer)+address);
                break;
            case NIFTI_TYPE_FLOAT128:
            case NIFTI_TYPE_COMPLEX128:
            case NIFTI_TYPE_COMPLEX256:
            case NIFTI_TYPE_COMPLEX64:
            default:
                return FALSE;
        }
    }
    
    rsScaleValues(timecourse, th, slope, inter);
    
    return TRUE;
}

/*
 * Takes in a buffer obtained from fsl and reads
 * out the values for all the voxels that are
 * specified by the parameter points for a
 * certain timepoint t.
 * The variables xh, yh, zh, and th correspond
 * to the dimensions of the nifti file.
 *
 */
BOOL rsExtractPointsFromBuffer(const FSLIO *fslio, double *data, const void *buffer, const float slope, const float inter, const Point3D* points, const unsigned long nPoints, const int t, const int xh, const int yh, const int zh, const int th) {
    
    for (unsigned long iPoint=0; iPoint<nPoints; iPoint=iPoint+1) {
        
        const Point3D p = points[iPoint];
        const long voxeLOffset = p.z * (xh * yh) + (p.y * xh) + p.x;
        const long timeOffset  = xh * yh * zh;
        const long address = t * timeOffset + voxeLOffset;
        
        switch (fslio->niftiptr->datatype) {
            case NIFTI_TYPE_UINT8:
                data[iPoint] = (double) *((THIS_UINT8 *)(buffer)+address);
                break;
            case NIFTI_TYPE_INT8:
                data[iPoint] = (double) *((THIS_INT8 *)(buffer)+address);
                break;
            case NIFTI_TYPE_UINT16:
                data[iPoint] = (double) *((THIS_UINT16 *)(buffer)+address);
                break;
            case NIFTI_TYPE_INT16:
                data[iPoint] = (double) *((THIS_INT16 *)(buffer)+address);
                break;
            case NIFTI_TYPE_UINT64:
                data[iPoint] = (double) *((THIS_UINT64 *)(buffer)+address);
                break;
            case NIFTI_TYPE_INT64:
                data[iPoint] = (double) *((THIS_INT64 *)(buffer)+address);
                break;
            case NIFTI_TYPE_UINT32:
                data[iPoint] = (double) *((THIS_UINT32 *)(buffer)+address);
                break;
            case NIFTI_TYPE_INT32:
                data[iPoint] = (double) *((THIS_INT32 *)(buffer)+address);
                break;
            case NIFTI_TYPE_FLOAT32:
                data[iPoint] = (double) *((THIS_FLOAT32 *)(buffer)+address);
                break;
            case NIFTI_TYPE_FLOAT64:
                data[iPoint] = (double) *((THIS_FLOAT64 *)(buffer)+address);
                break;
            case NIFTI_TYPE_FLOAT128:
            case NIFTI_TYPE_COMPLEX128:
            case NIFTI_TYPE_COMPLEX256:
            case NIFTI_TYPE_COMPLEX64:
            default:
                return FALSE;
        }
    }
    
    rsScaleValues(data, nPoints, slope, inter);
    
    return TRUE;
}

/*
 * Takes a timecourse and re-scales it before writing
 * it into a buffer that than later can be written out
 * into a nifti.
 * The point p represents the voxel coordinate for
 * which the timecourse will be written. The variables
 * xh, yh, zh, and th correspond to the dimensions of
 * the nifti file.
 */
BOOL rsWriteTimecourseToBuffer(const FSLIO *fslio, const double *timecourse, void *buffer, const float slope, const float inter, const Point3D p, const int xh, const int yh, const int zh, const int th) {
    
    THIS_UINT8   *THIS_UINT8_BUFFER   = (THIS_UINT8*)buffer;
    THIS_INT8    *THIS_INT8_BUFFER    = (THIS_INT8*)buffer;
    THIS_UINT16  *THIS_UINT16_BUFFER  = (THIS_UINT16*)buffer;
    THIS_INT16   *THIS_INT16_BUFFER   = (THIS_INT16*)buffer;
    THIS_UINT32  *THIS_UINT32_BUFFER  = (THIS_UINT32*)buffer;
    THIS_INT32   *THIS_INT32_BUFFER   = (THIS_INT32*)buffer;
    THIS_UINT64  *THIS_UINT64_BUFFER  = (THIS_UINT64*)buffer;
    THIS_INT64   *THIS_INT64_BUFFER   = (THIS_INT64*)buffer;
    THIS_FLOAT32 *THIS_FLOAT32_BUFFER = (THIS_FLOAT32*)buffer;
    THIS_FLOAT64 *THIS_FLOAT64_BUFFER = (THIS_FLOAT64*)buffer;
    
    for (int t=0; t<th; t=t+1) {
        
        long address = t * (xh * yh * zh) + p.z * (xh * yh) + (p.y * xh) + p.x;
        
        switch (fslio->niftiptr->datatype) {
            case NIFTI_TYPE_UINT8:
                THIS_UINT8_BUFFER[address]   = (THIS_UINT8)   ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_INT8:
                THIS_INT8_BUFFER[address]    = (THIS_INT8)    ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_UINT16:
                THIS_UINT16_BUFFER[address]  = (THIS_UINT16)  ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_INT16:
                THIS_INT16_BUFFER[address]   = (THIS_INT16)   ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_UINT32:
                THIS_UINT32_BUFFER[address]  = (THIS_UINT32)  ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_INT32:
                THIS_INT32_BUFFER[address]   = (THIS_INT32)   ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_UINT64:
                THIS_UINT64_BUFFER[address]  = (THIS_UINT64)  ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_INT64:
                THIS_INT64_BUFFER[address]   = (THIS_INT64)   ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_FLOAT32:
                THIS_FLOAT32_BUFFER[address] = (THIS_FLOAT32) ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_FLOAT64:
                THIS_FLOAT64_BUFFER[address] = (THIS_FLOAT64) ((timecourse[t] - inter) / slope);
                break;
            case NIFTI_TYPE_FLOAT128:
            case NIFTI_TYPE_COMPLEX128:
            case NIFTI_TYPE_COMPLEX256:
            case NIFTI_TYPE_COMPLEX64:
            default:
                return FALSE;
        }
    }
    
    return TRUE;
}

size_t rsWordLength(const int datatype) {
	switch (datatype) {
        case NIFTI_TYPE_UINT8:
            return sizeof(THIS_UINT8);
            break;
        case NIFTI_TYPE_INT8:
            return sizeof(THIS_INT8);
            break;
        case NIFTI_TYPE_UINT16:
            return sizeof(THIS_UINT16);
            break;
        case NIFTI_TYPE_INT16:
            return sizeof(THIS_INT16);
            break;
        case NIFTI_TYPE_UINT32:
            return sizeof(THIS_UINT32);
            break;
        case NIFTI_TYPE_INT32:
            return sizeof(THIS_INT32);
            break;
        case NIFTI_TYPE_UINT64:
            return sizeof(THIS_UINT64);
            break;
        case NIFTI_TYPE_INT64:
			return sizeof(THIS_INT64);
            break;
        case NIFTI_TYPE_FLOAT32:
            return sizeof(THIS_FLOAT32);
            break;
        case NIFTI_TYPE_FLOAT64:
			return sizeof(THIS_FLOAT64);
            break;
        case NIFTI_TYPE_FLOAT128:
        case NIFTI_TYPE_COMPLEX128:
        case NIFTI_TYPE_COMPLEX256:
        case NIFTI_TYPE_COMPLEX64:
        default:
            return 0;
    }
}

size_t rsVolumeOffset(const int xh, const int yh, const int zh, const int t) {
	return (size_t)t * ((size_t)xh * (size_t)yh * (size_t)zh);
}

size_t rsVolumeLength(const int xh, const int yh, const int zh) {
	return rsVolumeOffset(xh, yh, zh, 1);
}

/*
 * Takes in a buffer obtained from fsl and reads
 * out the values for one volume.
 *
 */
BOOL rsExtractVolumeFromBuffer(const FSLIO *fslio, double *data, const void *buffer, const float slope, const float inter, const int t, const int xh, const int yh, const int zh) {
	
	const int nifti_datatype  = fslio->niftiptr->datatype;
	const size_t wordlength   = rsWordLength(nifti_datatype);
	const size_t volumeOffset = rsVolumeOffset(xh, yh, zh, t);
	const size_t volumeLength = rsVolumeLength(xh, yh, zh);
	
	void *volBuffer = (void*)((char*)buffer + wordlength*volumeOffset);
	
	return convertBufferToScaledDouble(data, volBuffer, volumeLength, slope, inter, nifti_datatype);
}

/*
 * Takes in a buffer obtained from fsl and reads
 * out the values for one volume.
 *
 */
BOOL rsWriteVolumeToBuffer(const FSLIO *fslio, double *data, const void *buffer, const float slope, const float inter, const int t, const int xh, const int yh, const int zh) {
	
	const int nifti_datatype  = fslio->niftiptr->datatype;
	const size_t wordlength   = rsWordLength(nifti_datatype);
	const size_t volumeOffset = rsVolumeOffset(xh, yh, zh, t);
	
	void *volBuffer = (void*)((char*)buffer + wordlength*volumeOffset);
	
	return convertScaledDoubleToBuffer(nifti_datatype, volBuffer, data, slope, inter, xh, yh, zh, FALSE);
}


/*
 * This function will convert a given 3D matrix of scaled doubles into
 * a buffer using the supplied datatype.
 * In the conversion the data will be rescaled(inter and slope) as
 * required by FslWriteVolumes.
 * It is the counterpart to the method convertBufferToScaledDouble that
 * is included in FslIO.
 */
BOOL convertScaledDoubleToBuffer(int datatype, void *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    switch ( datatype ) {
        case NIFTI_TYPE_UINT8:
            convertScaledDoubleToBuffer_UINT8(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_INT8:
            convertScaledDoubleToBuffer_INT8(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_UINT16:
            convertScaledDoubleToBuffer_UINT16(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_INT16:
            convertScaledDoubleToBuffer_INT16(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_UINT64:
            convertScaledDoubleToBuffer_UINT64(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_INT64:
            convertScaledDoubleToBuffer_INT64(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_UINT32:
            convertScaledDoubleToBuffer_UINT32(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_INT32:
            convertScaledDoubleToBuffer_INT32(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_FLOAT32:
            convertScaledDoubleToBuffer_FLOAT32(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
        case NIFTI_TYPE_FLOAT64:
            convertScaledDoubleToBuffer_FLOAT64(outbuf, inbuf, slope, inter, xh, yh, zh, multidim);
            return TRUE;
            
        case NIFTI_TYPE_FLOAT128:
        case NIFTI_TYPE_COMPLEX128:
        case NIFTI_TYPE_COMPLEX256:
        case NIFTI_TYPE_COMPLEX64:
        default:
            return FALSE;
            
    }
    
    return TRUE;
}

BOOL rsResetBufferToValue(const int datatype, void *buffer, const float slope, const float inter, const int xh, const int yh, const int zh, const int th, const double value) {
    
    switch (datatype) {
        case NIFTI_TYPE_UINT8:
        case NIFTI_TYPE_INT8:
        case NIFTI_TYPE_UINT16:
        case NIFTI_TYPE_INT16:
        case NIFTI_TYPE_UINT32:
        case NIFTI_TYPE_INT32:
        case NIFTI_TYPE_UINT64:
        case NIFTI_TYPE_INT64:
        case NIFTI_TYPE_FLOAT32:
        case NIFTI_TYPE_FLOAT64:
            break;
        case NIFTI_TYPE_FLOAT128:
        case NIFTI_TYPE_COMPLEX128:
        case NIFTI_TYPE_COMPLEX256:
        case NIFTI_TYPE_COMPLEX64:
        default:
            return FALSE;
    }

    
    THIS_UINT8   *THIS_UINT8_BUFFER   = (THIS_UINT8*)buffer;
    THIS_INT8    *THIS_INT8_BUFFER    = (THIS_INT8*)buffer;
    THIS_UINT16  *THIS_UINT16_BUFFER  = (THIS_UINT16*)buffer;
    THIS_INT16   *THIS_INT16_BUFFER   = (THIS_INT16*)buffer;
    THIS_UINT32  *THIS_UINT32_BUFFER  = (THIS_UINT32*)buffer;
    THIS_INT32   *THIS_INT32_BUFFER   = (THIS_INT32*)buffer;
    THIS_UINT64  *THIS_UINT64_BUFFER  = (THIS_UINT64*)buffer;
    THIS_INT64   *THIS_INT64_BUFFER   = (THIS_INT64*)buffer;
    THIS_FLOAT32 *THIS_FLOAT32_BUFFER = (THIS_FLOAT32*)buffer;
    THIS_FLOAT64 *THIS_FLOAT64_BUFFER = (THIS_FLOAT64*)buffer;
    
    double scaled_value = (value - inter) / slope;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) shared(buffer,THIS_UINT8_BUFFER,THIS_INT8_BUFFER,THIS_UINT16_BUFFER,THIS_INT16_BUFFER,THIS_UINT32_BUFFER,THIS_INT32_BUFFER,THIS_UINT64_BUFFER,THIS_INT64_BUFFER,THIS_FLOAT32_BUFFER,THIS_FLOAT64_BUFFER)
    {
        #pragma omp for schedule(guided)
        for (int x=0; x<xh; x=x+1) {
            for (int y=0; y<yh; y=y+1) {
                for (int z=0; z<zh; z=z+1) {
                    for (int t=0; t<th; t=t+1) {
                        long address = t * (xh * yh * zh) + z * (xh * yh) + (y * xh) + x;
                        
                        switch (datatype) {
                            case NIFTI_TYPE_UINT8:
                                THIS_UINT8_BUFFER[address]   = (THIS_UINT8)   scaled_value;
                                break;
                            case NIFTI_TYPE_INT8:
                                THIS_INT8_BUFFER[address]    = (THIS_INT8)    scaled_value;
                                break;
                            case NIFTI_TYPE_UINT16:
                                THIS_UINT16_BUFFER[address]  = (THIS_UINT16)  scaled_value;
                                break;
                            case NIFTI_TYPE_INT16:
                                THIS_INT16_BUFFER[address]   = (THIS_INT16)   scaled_value;
                                break;
                            case NIFTI_TYPE_UINT32:
                                THIS_UINT32_BUFFER[address]  = (THIS_UINT32)  scaled_value;
                                break;
                            case NIFTI_TYPE_INT32:
                                THIS_INT32_BUFFER[address]   = (THIS_INT32)   scaled_value;
                                break;
                            case NIFTI_TYPE_UINT64:
                                THIS_UINT64_BUFFER[address]  = (THIS_UINT64)  scaled_value;
                                break;
                            case NIFTI_TYPE_INT64:
                                THIS_INT64_BUFFER[address]   = (THIS_INT64)   scaled_value;
                                break;
                            case NIFTI_TYPE_FLOAT32:
                                THIS_FLOAT32_BUFFER[address] = (THIS_FLOAT32) scaled_value;
                                break;
                            case NIFTI_TYPE_FLOAT64:
                                THIS_FLOAT64_BUFFER[address] = (THIS_FLOAT64) scaled_value;
                                break;
                        }
                    }
                }
            }
        }
    }
    
    return TRUE;
}

void convertScaledDoubleToBuffer_UINT8(THIS_UINT8 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;

    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_UINT8) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_INT8(THIS_INT8 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;

    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_INT8) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_UINT16(THIS_UINT16 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;

    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_UINT16) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_INT16(THIS_INT16 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;

    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_INT16) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_UINT64(THIS_UINT64 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;

    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_UINT64) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_INT64(THIS_INT64 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;

    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_INT64) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_UINT32(THIS_UINT32 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;
    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_UINT32) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_INT32(THIS_INT32 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;
    
    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_INT32) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_FLOAT32(THIS_FLOAT32 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;

    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_FLOAT32) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_FLOAT64(THIS_FLOAT64 *outbuf, double *inbuf, float slope, float inter, int xh, int yh, int zh, BOOL multidim) {
    int x, y, z;

    for (z=0; z<zh; z++) {
        for (y=0; y<yh; y++) {
            long inAddress = z * ((yh+multidim) * (xh+multidim)) + y * (xh+multidim);
            long outAddress = z * (yh * xh) + y * xh;
            for (x=0; x<xh; x++) {
                outbuf[outAddress + x] = (THIS_FLOAT64) ((inbuf[inAddress + x] - inter) / slope);
            }
        }
    }
}

char *rsLeftTrimString(char *s)
{
    while(isspace(*s)) s++;
    return s;
}

char *rsRightTrimString(char *s)
{
    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}

char *rsTrimString(char *s)
{
    return rsRightTrimString(rsLeftTrimString(s));
}

void rsWriteNiftiHeader(FSLIO *fslio, char* comment)
{
  short sform_code, qform_code;
  mat44 smat, qmat;
  /* writes header and opens img file ready for writing */
  if (fslio==NULL)  RSIOERR("FslWriteHeader: Null pointer passed for FSLIO");

  if (fslio->niftiptr!=NULL) {

	char* oldComment = "";
	BOOL commentExistsAlready = FALSE;
	nifti1_extension *ext;
		
	// check extensions
	nifti_image *nim = fslio->niftiptr;

	if( nim->num_ext > 0 && nim->ext_list != NULL ) {
		
		// find extension
		ext = nim->ext_list;

		int c;
		for ( c = 0; c < nim->num_ext; c++ ){
			if ( ext->ecode == NIFTI_ECODE_COMMENT && ext->edata != NULL ) { 
				break;
			}
		    ext++;
		}

		if ( c < nim->num_ext ) {
			commentExistsAlready = TRUE;
			
			// read extension
			size_t size = ext->esize;
			oldComment = malloc(sizeof(char)*(size+2));
			sprintf(&oldComment[0], "%s\n", ext->edata);
		}
	}
	
	// incorporate current date/time into the comment
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	char date[23];
	sprintf(&date[0], "# %04d/%02d/%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

	// incorporate rstools version number into the comment
	char *version = "\n# " RSTOOLS_VERSION_LABEL "\n";
	
	// merge old/new comment, version and date
	size_t datalength       = strlen(oldComment)+strlen(version)+strlen(comment)+strlen(date)+1;
	size_t datalength2      = (size_t)ceil((double)datalength / 16.0) * 16;
	char data[datalength2];
	
	sprintf(&data[0], "%s%s%s%s", oldComment, version, date, comment);
	
	if ( commentExistsAlready ) {
		ext->edata = data;
		ext->esize = datalength2;
	} else {
		nifti_add_extension(fslio->niftiptr, &data[0], datalength2, NIFTI_ECODE_COMMENT);
	}
	
    fslio->written_hdr = 1;
    if (znz_isnull(fslio->fileptr)) RSIOERR("FslWriteHeader: no file opened!");
    /* modify niftiptr for FSL-specific purposes */
    strcpy(fslio->niftiptr->descrip, RSTOOLS_VERSION_LABEL);

    /* set qform to equal sform if currently unset (or vice versa) */
    qform_code = FslGetRigidXform(fslio,&qmat);
    sform_code = FslGetStdXform(fslio,&smat);
    if ( (sform_code != NIFTI_XFORM_UNKNOWN) && 
	 (qform_code == NIFTI_XFORM_UNKNOWN) ) {
      FslSetRigidXform(fslio,sform_code,smat);
    }
    if ( (qform_code != NIFTI_XFORM_UNKNOWN) && 
	 (sform_code == NIFTI_XFORM_UNKNOWN) ) {
      FslSetStdXform(fslio,qform_code,qmat);
    }
    if (FslIsSingleFileType(FslGetFileType(fslio))) {
      /* write header info but don't close the file */
      nifti_image_write_hdr_img2(fslio->niftiptr,2,"wb",fslio->fileptr,NULL);
      /* set up pointer at end of iname_offset for single files only */
      FslSeekVolume(fslio,0);
    } else {
      /* open a new hdr file, write it and close it */
      nifti_image_write_hdr_img(fslio->niftiptr,0,"wb");
    }
  }

  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return;
}

char *rsMergeStringArray(int argc, char * argv[])
{
	size_t length = 1;
	size_t written = 0;
	for ( int i=0; i<argc; i=i+1 ) {
		length = length + strlen(argv[i]) + 1;
	}
	
	char *result     = malloc(sizeof(char)*length);
	char *cur_string = result;
	char *separator  = " ";
	
	for ( int i=0; i<argc; i=i+1 ) {
		size_t cur_length = strlen(argv[i]);

		cur_string      = &result[0]+written;
		written         = written + cur_length;
		
		sprintf(&cur_string[0], "%s", argv[i]);
		
		cur_string      = &result[0]+written;
		strcpy(cur_string, separator);
		written         = written + 1;
		cur_string      = &result[0]+written;
	}
	cur_string[0] = '\0';
	
	return result;
}

char *rsReadCommentFile(char *path) {
	
	/* Open file */
	FILE *f = fopen(path, "r");
    
    if (f == NULL) {
        fprintf(stderr, "Error: Comments could not be read(file couldn't be opened).\n");
        return NULL;
    }

    /* Determine file size */
	fseek(f, 0L, SEEK_END);
	size_t size = ftell(f);
	rewind(f);
    
    /* Read out file */
    char *content = malloc(sizeof(char)*(size+1));
		
	if( content == NULL || size != fread(content, 1, size, f) ) {
		fprintf(stderr, "Error: Comments could not be read into memory!\n");
		free(content);
        return NULL;
	}
	fclose(f);
	
	content[size] = '\0';
	
	return &content[0];
}