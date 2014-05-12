//
//  rsbandpass2.c
//  rstools
//
//  Created by André Hoffmann on 6/28/13.
//
//
#include "rsbandpass_common.h"

int show_help( )
{
    printf(
	   RSTOOLS_VERSION_LABEL "\n"
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
    fslioFiltered = FslOpen(p.saveFilteredPath, "wb");
        
    if (fslioFiltered == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", p.saveFilteredPath);
        return 1;
    }
    
    FslCloneHeader(fslioFiltered, p.fslio);
    FslSetDim(fslioFiltered, p.xDim, p.yDim, p.zDim, p.vDim);
    FslSetDimensionality(fslioFiltered, 4);
    FslSetDataType(fslioFiltered, p.dt);
	char *callString = rsMergeStringArray(argc, argv);
    rsWriteNiftiHeader(fslioFiltered, callString);
	free(callString);
    
    // Prepare empty timecourse
    double emptybuffer[p.vDim];
    
    for (int i=0; i<p.vDim; i=i+1){
        emptybuffer[i] = log(-1.0);
    }
    
    // Prepare buffer
	size_t buffsize = rsGetBufferSize(p.xDim, p.yDim, p.zDim, p.vDim, p.dt);
    void *buffer    = malloc(buffsize);
    
    if (buffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
    
    FslReadVolumes(p.fslio, buffer, p.vDim);
    short x,y,z, processedSlices = 0;
    double *signal;
    Point3D point;
    
    #pragma omp parallel num_threads(p.threads) private(z,y,x,signal,point) shared(p,emptybuffer,buffer,fslioFiltered,processedSlices)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<p.zDim; z=z+1) {
            for (y=0; y<p.yDim; y=y+1) {
                for (x=0; x<p.xDim; x=x+1) {
                    
                    point = MakePoint3D(x, y, z);
                    
                    /* If it's not in the mask skip it to improve the performance */
                    if (p.mask != NULL && p.mask[z][y][x] < 0.1) {
                    
                        /* set the value in the filtered data to NaN so that the nifti isn't empty */
                        rsWriteTimecourseToBuffer(fslioFiltered, emptybuffer, buffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.vDim);
                        continue;
                    }
                    
                    /* read out timecourse from buffer */
                    signal = malloc(p.vDim*sizeof(double));
                    rsExtractTimecourseFromBuffer(p.fslio, signal, buffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.vDim);
                    
                    /* apply filter */
                    rsFFTFilter(p.fftParams, signal);
                    
                    /* write out filtered data to buffer */
                    rsWriteTimecourseToBuffer(fslioFiltered, signal, buffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.vDim);
                    free(signal);
                }
            }
			
			/* show progress */
			if (p.verbose) {
            	#pragma omp atomic
            	processedSlices += 1;
            
            	if (processedSlices > 0 && processedSlices % (short)(p.zDim / 10) == 0) {
                	fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)p.zDim));
            	}
			}
        }
    }
    
	if ( p.verbose ) {
    	fprintf(stdout, "Write out result to: %s\n", p.saveFilteredPath);
	}
    
    FslWriteVolumes(fslioFiltered, buffer, p.vDim);
    FslClose(fslioFiltered);
    free(fslioFiltered);
    
    if ( p.maskpath != NULL ) {
        free(p.mask);
    }
    
    free(buffer);
    FslClose(p.fslio);
    free(p.fslio);
    
    rsFFTFilterFree(p.fftParams);
    
	return 0;
}