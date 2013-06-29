#include "rsbandpass_common.h"

int show_help( )
{
    printf(
       "rsbandpass: Given a 4D-Nifti and a frequency band this tool will apply\n"
       "            FFT filtering on it.\n"
       "\n"
    );
    
    rsBandpassPrintHelp();
    
   return 0;
}

int main(int argc, char * argv[])
{
    
	if( argc < 2 ) return show_help();
	
	struct rsBandpassParameters p = rsBandpassLoadParams(argc, argv);

    if ( ! p.parametersValid ) {
        fprintf(stderr, "Invalid parameters\n");
        return 1;
    }

    /* prepare filtered file */
    FSLIO *fslioFiltered;
   	void *filteredBuffer;
    fslioFiltered = FslOpen(p.saveFilteredPath, "wb");
        
    if (fslioFiltered == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", p.saveFilteredPath);
        return 1;
    }
    
    FslCloneHeader(fslioFiltered, p.fslio);
    FslSetDim(fslioFiltered, p.xDim, p.yDim, p.zDim, p.vDim);
    FslSetDimensionality(fslioFiltered, 4);
    FslSetDataType(fslioFiltered, p.pixtype);
    FslWriteHeader(fslioFiltered);
    
    size_t buffsize = (size_t)((size_t)p.vDim*(size_t)p.dt/(size_t)8);
    filteredBuffer = malloc(buffsize);
        
    // Prepare buffer
    buffsize = (size_t)((size_t)p.vDim*(size_t)p.dt/(size_t)8);
    void *buffer = malloc(buffsize);
    double signal[p.vDim];
    
    /* Prepare empty timecourse */
    void *emptybuffer = malloc(buffsize);
    
    double v[p.vDim];
    for (short t=0; t<p.vDim; t=t+1) {
        v[t] = log(-1.0);
    }
    convertScaledDoubleToBuffer(fslioFiltered->niftiptr->datatype, emptybuffer, v, p.slope, p.inter, p.vDim, 1, 1, FALSE);
        
    /* Iterate over all voxels that are to be filtered */
    for (short z=0; z<p.zDim; z=z+1) {
        fprintf(stdout, "Filtering slice Z%03hd/%03hd\n", z+1, p.zDim);
        for (short y=0; y<p.yDim; y=y+1) {
            for (short x=0; x<p.xDim; x=x+1) {
                
                /* If it's not in the mask skip it to improve the performance */
                if (p.mask != NULL && p.mask[z][y][x] < 0.1) {
                    
                    /* set the value in the filtered data to 0 so that the nifti isn't empty */
                    if ( fslioFiltered != NULL ) {
                        rsWriteTimeSeries(fslioFiltered, emptybuffer, x, y, z, p.vDim);
                    }
                    
                    continue;
                }
                
                /* read out timecourse */
                FslReadTimeSeries(p.fslio, buffer, x, y, z, p.vDim);
                convertBufferToScaledDouble(signal, buffer, (long)p.vDim, p.slope, p.inter, p.fslio->niftiptr->datatype);
                
                /* apply filter */
                rsFFTFilter(p.fftParams, signal);
                
                /* write out filtered data */
                convertScaledDoubleToBuffer(fslioFiltered->niftiptr->datatype, filteredBuffer, signal, p.slope, p.inter, p.vDim, 1, 1, FALSE);
                rsWriteTimeSeries(fslioFiltered, filteredBuffer, x, y, z, p.vDim);
            }
        }
    }
    
    FslClose(fslioFiltered);
    free(fslioFiltered);
    free(filteredBuffer);
    
    if ( p.maskpath != NULL ) {
        free(p.mask);
    }
    
    free(buffer);
    FslClose(p.fslio);
    free(p.fslio);
    free(emptybuffer);
    
    rsFFTFilterFree(p.fftParams);
    
	return 0;
}