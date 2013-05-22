#include "rsniftiutils.h"

Point3D MakePoint3D(unsigned int x, unsigned int y, unsigned int z)
{
    Point3D a;
    a.x = x;
    a.y = y;
    a.z = z;
    return a;
}

void convertScaledDoubleToBuffer(int datatype, void *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    switch ( datatype ) {
        case NIFTI_TYPE_UINT8:
            convertScaledDoubleToBuffer_UINT8(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_INT8:
            convertScaledDoubleToBuffer_INT8(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_UINT16:
            convertScaledDoubleToBuffer_UINT16(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_INT16:
            convertScaledDoubleToBuffer_INT16(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_UINT64:
            convertScaledDoubleToBuffer_UINT64(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_INT64:
            convertScaledDoubleToBuffer_INT64(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_UINT32:
            convertScaledDoubleToBuffer_UINT32(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_INT32:
            convertScaledDoubleToBuffer_INT32(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_FLOAT32:
            convertScaledDoubleToBuffer_FLOAT32(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
        case NIFTI_TYPE_FLOAT64:
            convertScaledDoubleToBuffer_FLOAT64(outbuf, inbuf, slope, inter, xh, yh, zh);
            break;
            
        case NIFTI_TYPE_FLOAT128:
        case NIFTI_TYPE_COMPLEX128:
        case NIFTI_TYPE_COMPLEX256:
        case NIFTI_TYPE_COMPLEX64:
        default:
            fprintf(stderr, "\nWarning, %s not supported yet.\n",nifti_datatype_string(datatype));
            
    }
}


void convertScaledDoubleToBuffer_UINT8(THIS_UINT8 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_UINT8) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_INT8(THIS_INT8 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_INT8) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_UINT16(THIS_UINT16 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_UINT16) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_INT16(THIS_INT16 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_INT16) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_UINT64(THIS_UINT64 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_UINT64) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_INT64(THIS_INT64 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_INT64) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_UINT32(THIS_UINT32 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_UINT32) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_INT32(THIS_INT32 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_INT32) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_FLOAT32(THIS_FLOAT32 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_FLOAT32) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

void convertScaledDoubleToBuffer_FLOAT64(THIS_FLOAT64 *outbuf, double ***inbuf, float slope, float inter, int xh, int yh, int zh) {
    int x, y, z;

    for (x=0; x<xh; x++) {
        for (y=0; y<yh; y++) {
            for (z=0; z<zh; z++) {
                outbuf[z * (yh * xh) + y * xh + x] = (THIS_FLOAT64) ((inbuf[z][y][x] - inter) / slope);
            }
        }
    }
}

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