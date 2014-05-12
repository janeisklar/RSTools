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
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"

#define RSTIMECOURSE_ALGORITHM_MEAN   1
#define RSTIMECOURSE_ALGORITHM_PCA    2
#define RSTIMECOURSE_ALGORITHM_TPCA   3
#define RSTIMECOURSE_ALGORITHM_CTP    4
#define RSTIMECOURSE_ALGORITHM_STDDEV 5

void rsWriteSpatialMap(char *file, FSLIO *reference, Point3D *points, gsl_matrix *maps);
void rsTimecourseTest();

int show_help( void )
{
   printf(	
	  RSTOOLS_VERSION_LABEL "\n"
      "rstimecourse2: Given a 4D-Nifti, this tool extracts the time course\n"
      "               for a single voxel or the meaned average of a region\n"
      "               specified by a binary mask\n"
      "\n"
   );
    
   printf(
      "basic usage:  rstimecourse [-m <mask> [-a <algorithm>] [-savemask <mask>]] [-p <X> <Y> <Z>] -input <volume>\n"
      "\n"
   );
    
   printf(
      "options:\n"
      "  -a[lgorithm] <algorithm>   : the algorithm used to aggregate the data within\n"
      "                               a ROI, e.g. mean, stddev, pca, tpca or ctp\n"
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
      "  -mask2 <mask>              : (use only with -a ctp) a second mask specifying the ROI\n"
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
	  "  -retainComponents <int>    : (use only with -a [t]pca or ctp) number of PCA/CTP components\n"
	  "                               that will be outputted. This should be a multiple of two when\n"
	  "                               running CTP as the components will be distributed equally over\n"
	  "                               both masks.\n"
   );

   printf(
	  "  -useStandardScores         : (use only with -a [t]pca) remove mean and set std. dev to 1\n"
	  "                               prior to running pca\n"
   );

   printf(
	  "  -spatialMap <volume>       : (use only with -a [t]pca or ctp) store spatial map that is\n"
	  "                               created using the PCA/CTP components\n"
   );


	printf(
	  "  -eigenvalues <file>        : (use only with -a [t]pca or ctp) write out all eigenvalues\n"
	  "                               to the file that is specified with this option\n"
	);

    printf(                        
      "  -threads <int>             : number of threads used for processing\n"
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
	char *mask2path = NULL;
    char *savemaskpath = NULL;
	char *spatialmappath = NULL;
	char *eigenvaluespath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	size_t pixtype;
	short dt;
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
		} else if ( ! strcmp(argv[ac], "-mask2") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -mask2\n");
				return 1;
			}
			mask2path = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strncmp(argv[ac], "-m", 2) ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-eigenvalues") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -eigenvalues\n");
				return 1;
			}
			eigenvaluespath = argv[ac];  /* no string copy, just pointer assignment */
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
			} else if ( ! strcmp(argv[ac], "ctp") ) {
				algorithm = RSTIMECOURSE_ALGORITHM_CTP;
			} else if ( ! strcmp(argv[ac], "stddev") ) {
				algorithm = RSTIMECOURSE_ALGORITHM_STDDEV;
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
		} else if ( ! strcmp(argv[ac], "-test") ) {
			rsTimecourseTest();
			return 0;
    	} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	rsSetThreadsNum(threads);
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified!\n");
		return 1;
	}
	
	if ( maskpath == NULL && (x<0||y<0||z<0) ) {
		fprintf(stderr, "Either a binary mask or a voxel coordinate must be specified!\n");
		return 1;
	}
	
	if ( mask2path == NULL && algorithm == RSTIMECOURSE_ALGORITHM_CTP ) {
		fprintf(stderr, "A second binary mask must be supplied to run CTP! (use -mask2)\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Input file: %s\n", inputpath);
        fprintf(stdout, "Mask file: %s\n", maskpath);
		if ( x >= 0 )
	    	fprintf(stdout, "Voxel: %d %d %d\n", x, y, z);
		if ( algorithm == RSTIMECOURSE_ALGORITHM_MEAN ) {
        	fprintf(stdout, "Algorithm: Mean\n");
		} else if ( algorithm == RSTIMECOURSE_ALGORITHM_STDDEV ) {
			fprintf(stdout, "Algorithm: Standard Deviation\n");
		} else if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA || algorithm == RSTIMECOURSE_ALGORITHM_TPCA || algorithm == RSTIMECOURSE_ALGORITHM_CTP ) {
			if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA ) {
				fprintf(stdout, "Algorithm: sPCA\n");
			} else if ( algorithm == RSTIMECOURSE_ALGORITHM_TPCA ) {
				fprintf(stdout, "Algorithm: tPCA\n");
			} else if ( algorithm == RSTIMECOURSE_ALGORITHM_CTP ) {
				fprintf(stdout, "Algorithm: CTP\n");
			}
			if ( spatialmappath != NULL ) {
				fprintf(stdout, "Spatial map file: %s\n", spatialmappath);
			}
			if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA || algorithm == RSTIMECOURSE_ALGORITHM_TPCA ) {
				fprintf(stdout, "Variance to be retained: %.2f\n", minVariance);
			}
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
	pixtype = FslGetDataType(fslio, &dt);
	
    if ( verbose ) {
        fprintf(stdout, "Dt: %d Pixtype: %ld\n", dt, pixtype);
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
        
        buffsize = rsGetBufferSize(1, 1, 1, vDim, dt);
        
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
		unsigned long nPointsMask1 = 0L;
        Point3D *maskPoints = ReadMask(maskpath, xDim, yDim, zDim, &nPointsMask1, savemaskpath, fslio, NULL);
        if ( maskPoints == NULL) {
            fprintf(stderr, "\nError: Mask invalid.\n");
            FslClose(fslio);
            return 1;
        }
		nPoints = nPoints + nPointsMask1;
		
		/* Load second mask (for CST) */
		unsigned long nPointsMask2 = 0L;
		Point3D *mask2Points;
		if ( algorithm == RSTIMECOURSE_ALGORITHM_CTP ) {
	        if ( mask2path == NULL) {
	            fprintf(stderr, "\nError: Mask2 invalid.\n");
	            FslClose(fslio);
	            return 1;
	        } else {
				mask2Points= ReadMask(mask2path, xDim, yDim, zDim, &nPointsMask2, NULL, fslio, NULL);
			
				// merge with points of first mask
				nPoints = nPoints + nPointsMask2;
				Point3D* allPoints = malloc(nPoints*sizeof(Point3D));
				for ( long i=0; i<nPoints; i=i+1 ) {
					if ( i < nPointsMask1 ) {
						allPoints[i] = maskPoints[i];
					} else {
						allPoints[i] = mask2Points[i-nPointsMask1];
					}
				}
			
				free(mask2Points);
				free(maskPoints);
				maskPoints = allPoints;
			}
		}
        
        double timecourse[vDim];
        
        /* Read in volume */
        buffsize = rsGetBufferSize(xDim, yDim, zDim, vDim, dt);
        buffer = malloc(buffsize);
        FslReadVolumes(fslio, buffer, vDim);
        
		/* Identify irrelevant voxels */
		if ( verbose ) {
			fprintf(stdout, "Removing voxels that are either NaN or have a StdDev of 0\n");
		}
		
		long n;
		long nNanPoints = 0;
		long nNanPointsMask1 = 0;
		long nNanPointsMask2 = 0;
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
					if (n<nPointsMask1) {
						nNanPointsMask1 = nNanPointsMask1 + 1;
					} else {
						nNanPointsMask2 = nNanPointsMask2 + 1;
					}
					break;
				}
	        }

			if ( isNanPoint[n] == FALSE && gsl_stats_sd(signalData, 1, vDim) < 0.000001 ) {
				isNanPoint[n] = TRUE;
				nNanPoints = nNanPoints + 1;
				if (n<nPointsMask1) {
					nNanPointsMask1 = nNanPointsMask1 + 1;
				} else {
					nNanPointsMask2 = nNanPointsMask2 + 1;
				}
			}

	        free(signalData);
		}
		
		if ( verbose ) {
			fprintf(stdout, "Removed %ld voxels, %ld remaining\n", nNanPoints, (nPoints - nNanPoints));
		}

		if ( algorithm == RSTIMECOURSE_ALGORITHM_MEAN || algorithm == RSTIMECOURSE_ALGORITHM_STDDEV ) {
        
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
					const double length = nPoints - nNanPoints;
                
	                /* Read out datapoints from the buffer */
	                double *pointValues = malloc(nPoints*sizeof(double));
	                rsExtractPointsFromBuffer(fslio, pointValues, buffer, slope, inter, maskPoints, nPoints, t, xDim, yDim, zDim, vDim);
                
	                /* Iterate over all points in the mask and compute mean */
	                for ( unsigned long i=0; i<nPoints; i=i+1) {
						if ( isNanPoint[i] ) {
							continue;
						}
	                    sum = sum + pointValues[i];
	                }
	                timecourse[t] = sum / length;
                
					/* If necessary iterate over them once more to compute the standard scores(std. dev) */
					if ( algorithm == RSTIMECOURSE_ALGORITHM_STDDEV ) {
						
						sum = 0.0;

		                for ( unsigned long i=0; i<nPoints; i=i+1) {
							if ( isNanPoint[i] ) {
								continue;
							}
							sum = sum + pow(pointValues[i]-timecourse[t], 2.0) / (length-1.0);
					    }

					    timecourse[t] = sqrt(sum);
					}

	                free(pointValues);
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
			long nEigenvalues = (nPoints - nNanPoints);
			
			if ( algorithm == RSTIMECOURSE_ALGORITHM_PCA ) {
				pcaResult = rsPCA(data, minVariance, nComponents, verbose);
				components = pcaResult.transformed;
				
				if ( spatialmappath != NULL) {
					rsWriteSpatialMap(spatialmappath, fslio, nonNanPoints, pcaResult.eigenvectors);
				}
			} else {
				nEigenvalues = vDim;
				pcaResult = rsTPCA(data, minVariance, nComponents, verbose);
				components = pcaResult.eigenvectors;
				
				if ( spatialmappath != NULL) {
					rsWriteSpatialMap(spatialmappath, fslio, nonNanPoints, pcaResult.transformed);
				}
			}
			gsl_matrix_free(data);
			free(nonNanPoints);
			
			// Write out eigenvalues
			if ( eigenvaluespath != NULL ) {
				FILE *fp;
				fp = fopen( eigenvaluespath , "w" );
				
				for ( long i=0; i<nEigenvalues; i=i+1 ) {
					fprintf(fp, "%.10f\n", gsl_vector_get(pcaResult.eigenvalues_all, i));
				}
				
				fclose(fp);
			}
			
			
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
					
		} else if ( algorithm == RSTIMECOURSE_ALGORITHM_CTP ) {
		
			if ( verbose ) {
				fprintf(stdout, "Computing common temporal patterns for the timecourses in both masks\n");
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
			gsl_matrix* data = gsl_matrix_alloc(vDim, (nPoints - nNanPoints));
			#pragma omp parallel num_threads(rsGetThreadsNum()) shared(fslio,buffer,slope,inter,xDim,yDim,zDim,vDim,nPoints,nNanPoints) private(n)
			{
				#pragma omp for schedule(guided)
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
						gsl_matrix_set(data, t, n, (signalData[t] - mean) / stdev);
					}
			
					free(signalData);
				}
			}
			free(buffer);
			
			// Run CTP
			gsl_matrix_view viewA = gsl_matrix_submatrix(data, 0,                            0, vDim, nPointsMask1-nNanPointsMask1);
			gsl_matrix_view viewB = gsl_matrix_submatrix(data, 0, nPointsMask1-nNanPointsMask1, vDim, nPointsMask2-nNanPointsMask2);
			gsl_matrix *dataA = &(viewA.matrix);
			gsl_matrix *dataB = &(viewB.matrix);
			//long nEigenvalues = (nPoints - nNanPoints);
			struct rsCTPResult ctpResult = rsCTP(dataA, dataB, nComponents/2, verbose);
			long nEigenvalues = ctpResult.eigenvalues_all->size;
			
			if ( spatialmappath != NULL) {
				//rsWriteSpatialMap(spatialmappath, fslio, nonNanPoints, pcaResult.eigenvectors);
			}

			gsl_matrix_free(data);
			free(nonNanPoints);
		
			// Write out eigenvalues
			if ( eigenvaluespath != NULL ) {
				FILE *fp;
				fp = fopen( eigenvaluespath , "w" );
			
				for ( long i=0; i<nEigenvalues; i=i+1 ) {
					fprintf(fp, "%.10f\n", gsl_vector_get(ctpResult.eigenvalues_all, i));
				}
			
				fclose(fp);
			}
		
			if ( verbose ) {
				fprintf(stdout, "Selected eigenvectors(%dx%d) of the CST anaylsis:\n", (int)ctpResult.eigenvectors->size1, (int)ctpResult.eigenvectors->size2);
			}
			for ( t = 0; t<vDim; t=t+1 ) {
				for ( int c = 0; c<ctpResult.eigenvectors->size1; c=c+1) {
	            	fprintf(stdout, "% 5.10f", gsl_matrix_get(ctpResult.eigenvectors, c, t));
					if ( c < (ctpResult.eigenvectors->size1-1) ) {
						fprintf(stdout, "\t");
					}
				}
				fprintf(stdout, "\n");
	        }

			rsCTPResultFree(ctpResult);
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
	short xDim, yDim, zDim, vDim, dt;
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
	FslGetDataType(reference, &dt);
	FslGetDim(reference, &xDim, &yDim, &zDim, &vDim);
	
	/* prepare file */
	FslCloneHeader(fslio, reference);
    FslSetDim(fslio, xDim, yDim, zDim, nMaps);
    FslSetDimensionality(fslio, 4);
    FslSetDataType(fslio, DT_FLOAT32);
    FslWriteHeader(fslio);
    
	/* prepare buffer */
    buffer = malloc(rsGetBufferSize(xDim, yDim, zDim, vDim, DT_FLOAT32));
	rsResetBufferToValue(reference->niftiptr->datatype, buffer, slope, inter, xDim, yDim, zDim, nMaps, sqrt(-1.0));
	
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

void rsTimecourseTest()
{
	/*
     >> a = rand(5);
     >> b = triu(a) + triu(a,1)';
     b =
     
     0.2760    0.4984    0.7513    0.9593    0.8407
     0.4984    0.9597    0.2551    0.5472    0.2543
     0.7513    0.2551    0.5060    0.1386    0.8143
     0.9593    0.5472    0.1386    0.1493    0.2435
     0.8407    0.2543    0.8143    0.2435    0.9293
    */
	
	double **m = d2matrix(4,4);
    
    // 1st row
    m[0][0] = 0.2760;
    m[0][1] = 0.4984; m[1][0] = 0.4984;
    m[0][2] = 0.7513; m[2][0] = 0.7513;
    m[0][3] = 0.9593; m[3][0] = 0.9593;
    m[0][4] = 0.8407; m[4][0] = 0.8407;
    
    // 2nd row
    m[1][1] = 0.9597;
    m[1][2] = 0.2551; m[2][1] = 0.2551;
    m[1][3] = 0.5472; m[3][1] = 0.5472;
    m[1][4] = 0.2543; m[4][1] = 0.2543;
    
    // 3rd row
    m[2][2] = 0.5060;
    m[2][3] = 0.1386; m[3][2] = 0.1386;
    m[2][4] = 0.8143; m[4][2] = 0.8143;
    
    // 4th row
    m[3][3] = 0.1493;
    m[3][4] = 0.2435; m[4][3] = 0.2435;
    
    // 5th row
    m[4][4] = 0.9293;

	gsl_matrix* A = gsl_matrix_alloc(5,5);

	for ( int r=0; r<5; r=r+1 ) {
		for ( int c=0; c<5; c=c+1 ) {
			gsl_matrix_set(A, r, c, m[r][c]);
		}
	}
	
	fprintf(stdout, "Original matrix:\n");	
	rs_gsl_matrix_fprintf(stdout, A, "%.4f");
	fprintf(stdout, "\n");
	
	long evChanged = rsMakePositiveDefiniteSymmetric(A);

	fprintf(stdout, "* changed %ld eigenvalues\n\n", evChanged);	
	
	fprintf(stdout, "Matrix after calling rsMakePositiveDefiniteSymmetric():\n");
	rs_gsl_matrix_fprintf(stdout, A, "%.4f");
	fprintf(stdout, "\n");
}