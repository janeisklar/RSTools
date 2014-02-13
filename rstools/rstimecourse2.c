//
//  rstimecourse2.c
//  rstools
//
//  Created by André Hoffmann on 9/11/13.
//
//

#include <stdio.h>
#include <strings.h>
#include <omp.h>

#include <nifti1.h>
#include <fslio.h>
#include "rsniftiutils.h"
#include "rsmathutils.h"

#define RSTIMECOURSE_ALGORITHM_MEAN 1
#define RSTIMECOURSE_ALGORITHM_PCA  2
#define RSTIMECOURSE_ALGORITHM_TPCA 3

void rsWriteSpatialMap(char *file, FSLIO *reference, Point3D *points, gsl_matrix *maps);

int show_help( void )
{
   printf(
      "rstimecourse2: Given a 4D-Nifti, this tool extracts the time course\n"
      "               for a single voxel or the meaned average of a region\n"
      "               specified by a binary mask\n"
      "\n"
   );
    
   printf(
      "basic usage:  rstimecourse2 [-m <mask> [-a <algorithm>] [-savemask <mask>]] [-p <X> <Y> <Z>] -input <volume>\n"
      "\n"
   );
    
   printf(
      "options:\n"
      "  -a[lgorithm] <algorithm>   : the algorithm used to aggregate the data within\n"
      "                               a ROI, e.g. mean, pca or tpca\n"
   );

   printf(
      "  -help                      : show this help\n"
   );
 
   printf(
      "  -input <volume>            : the volume from which the timecourse will be extracted\n"
   );
   
   printf(
      "  -mask <mask>               : a mask specifying the ROI\n"
   );
   
   printf(
      "  -p <X> <Y> <Z>             : speficies a voxel using nifti coordinates(0-based) from\n"
      "                               which the timecourse is to be extracted\n"
   );
    
   printf(
      "  -savemask <mask>           : optional path where the rescaled mask specified with -mask\n"
      "                               will be saved. The saved file with have the same dimensions\n"
      "                               as the input volume.\n"
   );

   printf(
	  "  -retainVariance <float>    : (use only with -a [t]pca) percentage of explained variance\n"
	  "                               that will be retained. keep in mind that a higher percentage\n"
	  "                               will result in more variables are to be returned.\n"
   );

   printf(
	  "  -retainComponents <int>    : (use only with -a [t]pca) number of PCA components that will be\n"
	  "                               outputted\n"
   );

   printf(
	  "  -useStandardScores         : (use only with -a [t]pca) remove mean and set std. dev to 1\n"
	  "                               prior to running pca\n"
   );

   printf(
	  "  -spatialMap <volume>       : (use only with -a pca) store spatial map that is created using\n"
	  "                               the PCA components\n"
   );
   
   printf(
      "  -v[erbose]                 : show debug information\n"
      "\n"
   );
    
   return 0;
}

int main(int argc, char * argv[])
{
	FSLIO *fslio;
	void *buffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *maskpath = NULL;
    char *savemaskpath = NULL;
	char *spatialmappath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
    float inter = 0.0, slope = 1.0;
    int threads = 1;
	float minVariance = 1.0;
	int nComponents = -1;
	short algorithm = RSTIMECOURSE_ALGORITHM_MEAN;
	BOOL useStandardScores = FALSE;
    BOOL verbose = FALSE;
	
	int ac;

	if( argc < 2 ) return show_help();   /* typing '-help' is sooo much work */

	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			return show_help();
		} else if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return 1;
			}
			inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-spatialMap") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -spatialMap\n");
				return 1;
			}
			spatialmappath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-s", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -savemask\n");
				return 1;
			}
			savemaskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-a", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -a[lgorithm]\n");
				return 1;
			}
			if ( ! strcmp(argv[ac], "mean") ) {
				algorithm = RSTIMECOURSE_ALGORITHM_MEAN;	
			} else if ( ! strcmp(argv[ac], "pca") ) {
				algorithm = RSTIMECOURSE_ALGORITHM_PCA;
			} else if ( ! strcmp(argv[ac], "tpca") ) {
				algorithm = RSTIMECOURSE_ALGORITHM_TPCA;
			} else {
				fprintf(stderr, "** the requested algorithm is not supported\n");
				return 1;
			}
		} else if ( ! strcmp(argv[ac], "-retainVariance") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -retainVariance\n");
           		return 1;
           	}
           	minVariance = atof(argv[ac]);
        } else if ( ! strcmp(argv[ac], "-retainComponents") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -retainComponents\n");
           		return 1;
           	}
           	nComponents = atoi(argv[ac]);
        } else if ( ! strcmp(argv[ac], "-useStandardScores") ) {
			useStandardScores = TRUE;
		} else if ( ! strncmp(argv[ac], "-v", 2) ) {
			verbose = TRUE;
		} else if ( ! strncmp(argv[ac], "-p", 2) ) {
			if( ac+3 >= argc ) {
				fprintf(stderr, "** missing argument for -p, 3 coordinates must be supplied!\n");
				return 1;
			}
			ac++;
			x = atoi(argv[ac]);
			ac++;
			y = atoi(argv[ac]);
			ac++;
			z = atoi(argv[ac]);
		} else if ( ! strcmp(argv[ac], "-threads") ) {
    		if( ++ac >= argc ) {
    			fprintf(stderr, "** missing argument for -threads\n");
    			return 1;
    		}
    		threads = atoi(argv[ac]);
    	} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified!\n");
		return 1;
	}
	
	if ( maskpath == NULL && (x<0||y<0||z<0) ) {
		fprintf(stderr, "Either a binary mask or a voxel coordinate must be specified!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file: %s\n", inputpath);
        fprintf(stdout, "Mask file: %s\n", maskpath);
		if ( x >= 0 )
	    	fprintf(stdout, "Voxel: %d %d %d\n", x, y, z);
		if ( algorithm == RSTIMECOURSE_ALGORITHM_MEAN ) {
        	fprintf(stdout, "Algorithm: Mean\n");
		} else if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA || algorithm == RSTIMECOURSE_ALGORITHM_TPCA ) {
			if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA ) {
				fprintf(stdout, "Algorithm: PCA\n");
				if ( spatialmappath != NULL ) {
					fprintf(stdout, "Spatial map file: %s\n", spatialmappath);
				}
			} else {
				fprintf(stdout, "Algorithm: tPCA\n");
			}
			fprintf(stdout, "Variance to be retained: %.2f\n", minVariance);
			if ( nComponents > 0 ) {
				fprintf(stdout, "Number of components to be retained: %d\n", nComponents);
			}
		}
    }
	
    fslio = FslOpen(inputpath, "rb");
    if (fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",inputpath);
        return 1;
    }

	/* determine dimensions */
	FslGetDim(fslio, &xDim, &yDim, &zDim, &vDim);
	
    if ( verbose ) {
        fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", xDim, yDim, zDim, vDim);
    }
    
    if (fslio->niftiptr->scl_slope != 0) {
        slope = fslio->niftiptr->scl_slope;
        inter = fslio->niftiptr->scl_inter;
    }
	
	/* determine datatype and initalize buffer */
	dt = FslGetDataType(fslio, &pixtype);
	
    if ( verbose ) {
        fprintf(stdout, "Dt: %ld Pixtype: %d\n", dt, pixtype);
    }
	
    if ( x >= 0 ) {
        /*** Extract timecourse for a single point ***/
        
        /* check inputs */
        if ( (x<0) || (x>=fslio->niftiptr->nx) ) {
            fprintf(stderr, "\nError: x index (%d) out of range [0..%d]\n",x,fslio->niftiptr->nx-1);
            FslClose(fslio);
            return 1;
        }
        if ( (y<0) || (y>=fslio->niftiptr->ny) ) {
            fprintf(stderr, "\nError: y index (%d) out of range [0..%d]\n",y,fslio->niftiptr->ny-1);
            FslClose(fslio);
            return 1;
        }
        if ( (z<0) || (z>=fslio->niftiptr->nz) ) {
            fprintf(stderr, "\nError: z index (%d) out of range [0..%d]\n",z,fslio->niftiptr->nz-1);
            FslClose(fslio);
            return 1;
        }
        
        buffsize = vDim*dt/8;
        
        if ( verbose ) {
            fprintf(stdout, "Buffsize: %ld\n", buffsize);
        }
        
        buffer = malloc(buffsize);
        
        /* read out timecourse */
        FslReadTimeSeries(fslio, buffer, x, y, z, vDim);
        
        for (t=0; t<vDim; t++) {
            void *newV = (char*)buffer + (t*dt/8);
            
            double v = 0.0;
            int ret = convertBufferToScaledDouble(&v, newV, 1, slope, inter, fslio->niftiptr->datatype);
         
            if ( verbose ) {
                fprintf(stdout, "%d: %.10f(%d)\n", t, v, ret);
            } else {
                fprintf(stdout, "%.10f\n", v);
            }
        }
        
        /* clear buffer */
        free(buffer);
    } else {
        /*** Extract timecourse for a mask ***/
        unsigned long nPoints = 0L;
        Point3D *maskPoints = ReadMask(maskpath, xDim, yDim, zDim, &nPoints, savemaskpath, fslio, NULL);
        if ( maskPoints == NULL) {
            fprintf(stderr, "\nError: Mask invalid.\n");
            FslClose(fslio);
            return 1;
        }
        
        double timecourse[vDim];
        
        /* Read in volume */
        buffsize = (size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)vDim*(size_t)dt/(size_t)8;
        buffer = malloc(buffsize);
        FslReadVolumes(fslio, buffer, vDim);
        
		/* Identify irrelevant voxels */
		if ( verbose ) {
			fprintf(stdout, "Removing voxels that are either NaN or have a StdDev of 0\n");
		}
		
		long n;
		long nNanPoints = 0;
		int isNanPoint[nPoints];
		
		for (n=0; n<nPoints; n=n+1) {
	        const Point3D point = maskPoints[n];
			isNanPoint[n] = FALSE;

	        // load signal
	        double *signalData = malloc(sizeof(double) * vDim);
	        rsExtractTimecourseFromBuffer(fslio, signalData, buffer, slope, inter, point, xDim, yDim, zDim, vDim);

	        for ( int t=0; t<vDim; t=t+1 ) {
				if ( signalData[t] != signalData[t] ) {
					isNanPoint[n] = TRUE;
					nNanPoints = nNanPoints + 1;
					break;
				}
	        }

			if ( gsl_stats_sd(signalData, 1, vDim) < 0.000001 ) {
				isNanPoint[n] = TRUE;
				nNanPoints = nNanPoints + 1;
			}

	        free(signalData);
		}
		
		if ( verbose ) {
			fprintf(stdout, "Removed %ld voxels, %ld remaining\n", nNanPoints, (nPoints - nNanPoints));
		}

		if ( algorithm == RSTIMECOURSE_ALGORITHM_MEAN ) {
        
	        int t;

			if ( verbose ) {
				fprintf(stdout, "Coputing mean ROI timecourse\n");
			}

	        /* Iterate over all timepoints */
	        #pragma omp parallel num_threads(threads) private(t) shared(timecourse,nPoints)
	        {
	            #pragma omp for schedule(guided)
	            for ( t = 0; t<vDim; t=t+1 ) {
            
	                double sum = 0;
                
	                /* Read out datapoints from the buffer */
	                double *pointValues = malloc(nPoints*sizeof(double));
	                rsExtractPointsFromBuffer(fslio, pointValues, buffer, slope, inter, maskPoints, nPoints, t, xDim, yDim, zDim, vDim);
                
	                /* Iterate over all points in the mask */
	                for ( unsigned long i=0; i<nPoints; i=i+1) {
						if ( isNanPoint[i] ) {
							continue;
						}
	                    sum = sum + pointValues[i];
	                }
                
	                free(pointValues);
                
	                /* Create average */
	                timecourse[t] = sum / (nPoints - nNanPoints);
	            }
	        }
        
        
	        for ( t = 0; t<vDim; t=t+1 ) {
	            fprintf(stdout, "%.10f\n", timecourse[t]);
	        }

		} else if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA || algorithm == RSTIMECOURSE_ALGORITHM_TPCA ) {
			
			if ( verbose ) {
				fprintf(stdout, "Computing PCA components for timecourses in the mask\n");
			}
			
			// Prepare non-NAN points
			Point3D *nonNanPoints = malloc((nPoints - nNanPoints)*sizeof(Point3D));
			n = 0;
			for ( unsigned long i=0; i<nPoints; i=i+1) {
				if ( isNanPoint[i] ) {
					continue;
				}
				
				nonNanPoints[n] = maskPoints[i];				
				n = n+1;
            }
			
			// Prepare data matrix
			gsl_matrix* data = gsl_matrix_alloc((nPoints - nNanPoints), vDim);
			for ( n=0; n<(nPoints - nNanPoints); n=n+1) {
				
				double *signalData = malloc(sizeof(double) * vDim);
		        rsExtractTimecourseFromBuffer(fslio, signalData, buffer, slope, inter, nonNanPoints[n], xDim, yDim, zDim, vDim);
				
				double mean  = 0; 
				double stdev = 1;
				
				if ( useStandardScores ) {
					mean  = gsl_stats_mean(signalData, 1, vDim);
					stdev = gsl_stats_sd(signalData, 1, vDim);
				}
				
				for ( int t = 0; t < vDim; t=t+1 ) {
					gsl_matrix_set(data, n, t, (signalData[t] - mean) / stdev);
				}
				
				free(signalData);
			}
			
			// Run PCA
			gsl_matrix* components;
			struct rsPCAResult pcaResult;
			
			if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA ) {
				pcaResult = rsPCA(data, minVariance, nComponents, verbose);
				components = pcaResult.transformed;
				
				if ( spatialmappath != NULL) {
					rsWriteSpatialMap(spatialmappath, fslio, nonNanPoints, pcaResult.eigenvectors);
				}
			} else {
				pcaResult = rsTPCA(data, minVariance, nComponents, verbose);
				components = pcaResult.eigenvectors;
				
				if ( spatialmappath != NULL) {
					rsWriteSpatialMap(spatialmappath, fslio, nonNanPoints, pcaResult.transformed);
				}
			}
			gsl_matrix_free(data);
			free(nonNanPoints);
			
			if ( verbose ) {
				if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA ) {
					fprintf(stdout, "Reduced data matrix(%dx%d) after spatial PCA:\n", (int)components->size1, (int)components->size2);
				} else {
					fprintf(stdout, "Selected tPCA eigenvectors(%dx%d)\n", (int)components->size1, (int)components->size2);
				}
			}
			for ( t = 0; t<vDim; t=t+1 ) {
				for ( int c = 0; c<components->size1; c=c+1) {
	            	fprintf(stdout, "% 5.10f", gsl_matrix_get(components, c, t));
					if ( c < (components->size1-1) ) {
						fprintf(stdout, "\t");
					}
				}
				fprintf(stdout, "\n");
	        }

			rsPCAResultFree(pcaResult);
		}
		
		free(maskPoints);
    }
    
    FslClose(fslio);
    free(fslio);
    
	return 0;
}

void rsWriteSpatialMap(char *file, FSLIO *reference, Point3D *points, gsl_matrix *maps)
{
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim, pixType;
	size_t buffsize;
    float inter = 0.0, slope = 1.0;
	void *buffer;
	short nMaps  = maps->size1;
	long nPoints = maps->size2;
	
	FSLIO *fslio = FslOpen(file, "wb");
    if (fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",file);
    }

    if (reference->niftiptr->scl_slope != 0) {
        slope = reference->niftiptr->scl_slope;
        inter = reference->niftiptr->scl_inter;
    }
	
	/* determine dimensions */
	FslGetDataType(reference, &pixType);
	FslGetDim(reference, &xDim, &yDim, &zDim, &vDim);
	
	/* prepare file */
	FslCloneHeader(fslio, reference);
    FslSetDim(fslio, xDim, yDim, zDim, nMaps);
    FslSetDimensionality(fslio, 4);
    FslSetDataType(fslio, DT_FLOAT32);
    FslWriteHeader(fslio);
    
	/* prepare buffer */
    buffer = malloc((size_t)xDim*(size_t)yDim*(size_t)zDim*(size_t)vDim*(size_t)reference->niftiptr->nbyper);
	rsResetBufferToValue(reference->niftiptr->datatype, buffer, slope, inter, xDim, yDim, zDim, nMaps, 1, sqrt(-1.0));
	
	/* write spatial maps to buffer */
	for (unsigned long p = 0L; p<nPoints; p=p+1L) {
		double data[nMaps];
		
		for ( short i=0; i<nMaps; i=i+1 ) {
			data[i] = gsl_matrix_get(maps, i, p);
		}
        
        rsWriteTimecourseToBuffer(fslio, &data[0], buffer, slope, inter, points[p], xDim, yDim, zDim, nMaps);
    }
    
	/* write to file */
    FslWriteVolumes(fslio, buffer, nMaps);
    FslClose(fslio);
    free(fslio);
}