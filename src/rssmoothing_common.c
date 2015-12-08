#include "rssmoothing_common.h"
#include "maths/utils.h"
#include "utils/rsio.h"
#include <math.h>

double ***rsCreateGaussianKernel(double sigma, short *xdimkernel, short *ydimkernel, short *zdimkernel, double xvoxsize, double yvoxsize, double zvoxsize);
void rsConvolveWithKernel(double ***result, double ***input, double ***kernel, short xdim, short ydim, short zdim, short xdimKernel, short ydimKernel, short zdimKernel);

void rsSmoothingInit(rsSmoothingParameters *p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->outputpath,
        RSIO_LASTFILE
    });
    
    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }
        
    /* open input file */
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }
    
    rsSetThreadsNum(p->threads);
    
    if ( p->kernelSizeFWHM < 0.00001 ) {
        fprintf(stderr, "\nError: The FWHM kernel size must be a positive floating point number! (got %.2f)\n", p->kernelSizeFWHM); 
        return;
    }
	
	p->kernelSigma = p->kernelSizeFWHM / (2.0*sqrt(2.0*log(2.0)));
	
    /* output the most important parameters to the user */
    if ( p->verbose ) {
        fprintf(stdout, "Input file:  %s\n", p->inputpath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
        fprintf(stdout, "FWHM kernel size:  %.5f mm\n", p->kernelSizeFWHM);
        fprintf(stdout, "Sigma: %.5f mm\n", p->kernelSigma);
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
    }
	
    p->parametersValid = TRUE;
}

void rsSmoothingRun(rsSmoothingParameters *p)
{
    p->parametersValid = FALSE;

	// get voxel spacing
	float xSpacing, ySpacing, zSpacing, tr;
	FslGetVoxDim(p->input->fslio, &xSpacing, &ySpacing, &zSpacing, &tr);

	// create gaussian kernel
	short kerneldim[3];
	double ***kernel = rsCreateGaussianKernel(p->kernelSigma, &kerneldim[0], &kerneldim[1], &kerneldim[2], xSpacing, ySpacing, zSpacing);
	double ***xKernel = d3matrix(0,0,kerneldim[0]-1);
	double ***yKernel = d3matrix(0,kerneldim[1]-1,0);
	double ***zKernel = d3matrix(kerneldim[2]-1,0,0);
	const short xMidKernel = (kerneldim[0]-1)/2,
	            yMidKernel = (kerneldim[1]-1)/2,
                zMidKernel = (kerneldim[2]-1)/2;
	for (short n=0; n<kerneldim[0]; n++)
		xKernel[0][0][n] = kernel[zMidKernel][yMidKernel][n];
	for (short n=0; n<kerneldim[1]; n++)
		yKernel[0][n][0] = kernel[zMidKernel][n][xMidKernel];
	for (short n=0; n<kerneldim[2]; n++)
		zKernel[n][0][0] = kernel[n][yMidKernel][xMidKernel];

    // create output file
    p->output = rsCloneNiftiFile(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, p->input->vDim);
    
    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the smoothed output (%s) could not be created.\n", p->outputpath);
        return;
    }
    
	rsWriteNiftiHeader(p->output->fslio, p->callString);

    // prepare the output file's content
    int t, processedVolumes = 0;
	double ***data;
	double ***tmp;
    omp_lock_t updateProgressLock;
    omp_init_lock(&updateProgressLock);
	
	#pragma omp parallel num_threads(rsGetThreadsNum()) private(t,data,tmp) shared(processedVolumes)
	{
	    #pragma omp for schedule(guided, 1)
	    for ( t=0; t<p->input->vDim; t++ ) {
	        // extract a single volume for timepoint t from the buffer
	        data = d3matrix(p->input->zDim-1, p->input->yDim-1, p->input->xDim-1);
			tmp  = d3matrix(p->input->zDim-1, p->input->yDim-1, p->input->xDim-1);
	        rsExtractVolumeFromRSNiftiFileBuffer(p->input, data[0][0], t);

			rsConvolveWithKernel(tmp, data, xKernel, p->input->xDim, p->input->yDim, p->input->zDim, kerneldim[0], 1, 1);
			rsConvolveWithKernel(data, tmp, yKernel, p->input->xDim, p->input->yDim, p->input->zDim, 1, kerneldim[1], 1);
			rsConvolveWithKernel(tmp, data, zKernel, p->input->xDim, p->input->yDim, p->input->zDim, 1, 1, kerneldim[2]);

	        // write back to the buffer
	        rsWriteVolumeToBuffer(p->output->dt, tmp[0][0], p->output->data, p->output->slope, p->output->inter, t, p->output->xDim, p->output->yDim, p->output->zDim);

	        // free up memory
	        free(data[0][0]); free(data[0]); free(data);
	        free(tmp[0][0]); free(tmp[0]); free(tmp);
			
            // show progress
            if (p->progressCallback != NULL) {
                omp_set_lock(&updateProgressLock);
                rsReportProgressEvent *event = (rsReportProgressEvent*)rsMalloc(sizeof(rsReportProgressEvent));
                event->run = processedVolumes;
                processedVolumes += 1;
                event->percentage = (double)processedVolumes*100.0 / (double)p->input->vDim;
                rsReportProgressCallback_t cb = p->progressCallback->cb;
                void *data = p->progressCallback->data;
                cb(event, data);
                rsFree(event);
                omp_unset_lock(&updateProgressLock);
            } else if (p->verbose) {
                #pragma omp atomic
                processedVolumes += 1;
            
                if (processedVolumes > 0 && processedVolumes % (short)(p->input->vDim >= 10 ? p->input->vDim / 10 : p->input->vDim) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedVolumes*100.0 / (float)p->input->vDim));
                }
            }
	    }
	}

    // write smoothed file
    FslWriteVolumes(p->output->fslio, p->output->data, p->input->vDim);

    p->parametersValid = TRUE;
}

void rsSmoothingDestroy(rsSmoothingParameters *p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->output != NULL ) {
        p->output->data = NULL;
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }

    rsSmoothingFreeParams(p);
}

void rsConvolveWithKernel(double ***result, double ***input, double ***kernel, short xdim, short ydim, short zdim, short xdimKernel, short ydimKernel, short zdimKernel) {
    const short midx=xdimKernel/2, 
	            midy=ydimKernel/2, 
				midz=zdimKernel/2;
    for (short z=0; z<zdim; z++) {
		for (short y=0; y<ydim; y++) {
			for (short x=0; x<xdim; x++) {
				size_t nValues = 0;
				double val=0.0, 
				       norm= 0.0;
				const short x3=x-midx, 
				            y3=y-midy, 
					        z3=z-midz;
				for (short mz=0; mz<zdimKernel; mz++) {
					for (short my=0; my<ydimKernel; my++) {
						for (short mx=0; mx<xdimKernel; mx++) {
							const short x2=x3+mx, 
							            y2=y3+my,
								        z2=z3+mz;
							if (x2>=0 && x2<xdim && y2>=0 && y2<ydim && z2>=0 && z2<zdim) {
								if (isnan(input[z2][y2][x2]) || isinf(input[z2][y2][x2])) {
									continue;
								}
								val+=input[z2][y2][x2] * kernel[mz][my][mx];
								norm+=kernel[mz][my][mx];
								nValues += 1;
	    					}
	  				  	}
					}
				}

				if (nValues > 0) {
					result[z][y][x] = fabs(norm) > 1e-12 ? val / norm : val;
				} else {
					// fill with NaN if we excluded all values in the kernel
					result[z][y][x] = log(-1.0);
				}
			}
		}
    }
}

double ***rsCreateGaussianKernel(double sigma, short *xdimkernel, short *ydimkernel, short *zdimkernel, double xvoxsize, double yvoxsize, double zvoxsize) {
	const double cutoff = 4.0,
                 norm   = 2*sigma*sigma,
				 dx2    = xvoxsize*xvoxsize,
				 dy2    = yvoxsize*yvoxsize,
				 dz2    = zvoxsize*zvoxsize;
	const short sx = ((short) ceil(sigma*cutoff/xvoxsize)),
	            sy = ((short) ceil(sigma*cutoff/yvoxsize)),
	            sz = ((short) ceil(sigma*cutoff/zvoxsize));
	*xdimkernel = 2*sx + 1;
	*ydimkernel = 2*sy + 1;
	*zdimkernel = 2*sz + 1;
	double ***kernel = d3matrix(*zdimkernel-1, *ydimkernel-1, *xdimkernel-1);
	
	for (short z=0; z<*zdimkernel; z++) {
		for (short y=0; y<*ydimkernel; y++) {
	    	for (short x=0; x<*xdimkernel; x++) {
				const short x2=x-sx,
				            y2=y-sy,
					        z2=z-sz;
				kernel[z][y][x]=exp(-(x2*x2*dx2+y2*y2*dy2+z2*z2*dz2)/norm);				
	    	}
	  	}
	}
	
	return kernel;
}


