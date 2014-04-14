//
//  rscorrelation2.c
//  rstools
//
//  Created by André Hoffmann on 7/9/13.
//
//

#include <omp.h>
#include "rscorrelation_common.h"
#include <gsl/gsl_cdf.h>

int show_help( void )
{
    rsCorrelationPrintHelp();
    
    printf(
       "   -tClusterSize <int>  : number of timepoints that will be clustered\n"
       "                          into one correlation map\n"
    );

    printf(
       "   -nClusters <int>     : number of timepoints clusters that will be\n"
       "                          computed and stored\n"
    );

    return 0;
}

void rsVaryingLengthWriteCorrelationFile(struct rsCorrelationParameters* p, void *correlationBuffer);
double* rsVaryingLengthZCorrelation(const double* X, const double* Y, const size_t length, const unsigned int nClusters, const unsigned int tClusteringSize);

int main(int argc, char * argv[]) {
	
	if( argc < 2 ) return show_help();
	
    struct rsCorrelationParameters p = rsCorrelationLoadParams(argc, argv);
    
    // Load file
    // Prepare buffer
    size_t buffsize = (size_t)p.xDim*(size_t)p.yDim*(size_t)p.zDim*(size_t)p.vDim*(size_t)p.dt/(size_t)8;
	void  *buffer;
	buffer = malloc(buffsize);
    
    if (buffer == NULL) {
        fprintf(stdout, "Not enough free memory :-(\n");
        return 1;
    }
	
	FslClose(p.fslioCorrelation);
	p.fslioCorrelation = FslOpen(p.outputpath, "wb");
	FslSetDim(p.fslioCorrelation, p.xDim, p.yDim, p.zDim, p.nClusters);
	FslSetIntent(fslioOutput, NIFTI_INTENT_ZSCORE, 0, 0, 0);
    rsWriteNiftiHeader(p.fslioCorrelation, p.comment);
    
    FslReadVolumes(p.fslio, buffer, p.vDim);
    
    // Prepare buffer
    short x,y,z, processedSlices = 0;
    double *timecourse;
	buffsize = (size_t)(rsVolumeLength(p.xDim, p.yDim, p.zDim)*(size_t)p.nClusters*(size_t)p.dt/(size_t)8);
    void *correlationBuffer;
    correlationBuffer = malloc(buffsize);
    
    /* Iterate over all voxels for which the correlation is to be computed */
    #pragma omp parallel num_threads(p.threads) private(y,x,timecourse) shared(p,processedSlices)
    {
        #pragma omp for schedule(guided)
        for (short z=0; z<p.zDim; z=z+1) {
            for (short y=0; y<p.yDim; y=y+1) {
                for (short x=0; x<p.xDim; x=x+1) {
                    
					Point3D point = MakePoint3D(x, y, z);

                    /* If it's not in the mask skip it to improve the performance */
                    if (p.mask != NULL && p.mask[z][y][x] < 0.1) {
                        
                        /* set the value in the correlation file to NaN so that it is skipped in later processing steps */
						double *correlations = malloc(sizeof(double)*p.nClusters);
						for (short c=0; c<p.nClusters; c=c+1) {
                        	correlations[c] = log(-1.0);
						}
						
                        rsWriteTimecourseToBuffer(p.fslioCorrelation, &correlations[0], correlationBuffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.nClusters);
						free(correlations);
						
                        continue;
                    }
                    
                    /* read out timecourse */
                    timecourse = malloc(p.vDim*sizeof(double));
                    rsExtractTimecourseFromBuffer(p.fslio, &timecourse[0], buffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.vDim);
                    
                    /* compute correlation */
					double *voxelCorrelations;
					voxelCorrelations = rsVaryingLengthZCorrelation(&timecourse[0], p.regressor, p.vDim, p.nClusters, p.timeClusteringSize);
					
					/* save correlations */
					rsWriteTimecourseToBuffer(p.fslioCorrelation, &voxelCorrelations[0], correlationBuffer, p.slope, p.inter, point, p.xDim, p.yDim, p.zDim, p.nClusters);

                    free(voxelCorrelations);
                    free(timecourse);
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
    
    /* Write correlation file */
    rsVaryingLengthWriteCorrelationFile(&p, correlationBuffer);
	free(correlationBuffer);

    /* Close everything */
    rsCorrelationFree(&p);
    free(buffer);
    
	return 0;
}

void rsVaryingLengthWriteCorrelationFile(struct rsCorrelationParameters* p, void *correlationBuffer) {
    
    /* Write correlation file */
    if ((*p).verbose) fprintf(stdout, "Writing correlation file\n");
    
    FslWriteVolumes((*p).fslioCorrelation, correlationBuffer, (*p).nClusters);
}

double* rsVaryingLengthZCorrelation(const double* X, const double* Y, const size_t length, const unsigned int nClusters, const unsigned int tClusteringSize) {
	
	// define some constants
    const long double epsilon = 0.0000000001;
    const long double half     = 0.5;
    const long double one      = 1.0;
    const long double minusone = -1.0;

	// define some working variables
	double *correlations;
	correlations = malloc(sizeof(double)*nClusters);
	double sampledX[length];
	double sampledY[length];
	size_t indices[length];

	// prepare array with indices that will be randomly drawn from
	for (size_t i = 0; i < length; i=i+1) {
		indices[i] = i;
	}

	// shuffle and draw samples
	gsl_ran_shuffle(rsGetRandomNumberGenerator(), indices, length, sizeof(size_t));

	//draw samples
	for (size_t i = 0; i < length; i=i+1) {
		const size_t index = indices[i];
		sampledX[i] = X[index];
		sampledY[i] = Y[index];
	}

	// compute correlations
	short N           = 1;
    short Nm1         = 0;
	short c           = 0;
	long double norm  = 0.0;
    long double sumX  = 0.0;
    long double sumY  = 0.0;
	long double meanX = 0.0;
	long double meanY = 0.0;
    long double tmpX  = 0.0;
    long double tmpY  = 0.0;
	long double tmpR  = 0.0;
    long double sX    = 0.0;
	long double sY    = 0.0;
    long double r     = 0.0;
	long double diffX = 0.0;
	long double diffY = 0.0;
    
	// increase the number of timepoints and update correlation
    for (size_t i=0; i<length; i=i+1) {
		// update sum and mean
        sumX  = sumX + (long double)sampledX[i]; meanX = sumX / (long double)N; diffX = sampledX[i]-meanX;
        sumY  = sumY + (long double)sampledY[i]; meanY = sumY / (long double)N; diffY = sampledY[i]-meanY;
		
		// update standard scores(std. dev)
		tmpX = tmpX + powl(diffX, 2.0);
		tmpY = tmpY + powl(diffY, 2.0);
		
	    // update correlation coefficient
        tmpR = tmpR + (diffX * diffY);		
		
		// end of next cluster reached?
		if ( (((i+1) % tClusteringSize) == 0 && i > 0) || (i+1) >= length ) {

			// normalize correlation coefficient
			sX   = sqrtl(tmpX / (long double)Nm1);
			sY   = sqrtl(tmpY / (long double)Nm1);
		    norm = fmaxl(sX, epsilon) * fmaxl(sY, epsilon) * (long double)Nm1;
			r = tmpR / norm;
			
			// fisher's r-to-z transformation
		    r = fminl(r, one-epsilon);
		    r = fmaxl(r, minusone+epsilon);
			r = half * logl( (one + r) / (one - r) );
			
			// store correlation
			correlations[c] = r;
			c = c+1;
			
			// max amount of clusters reached
			if ( c >= nClusters ) {
				break;
			}
		}
		
		// update length of timepoints
		Nm1 = N;
		N   = N + 1;
    }

	return &correlations[0];
}