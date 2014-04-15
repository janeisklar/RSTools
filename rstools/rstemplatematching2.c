//
//  rstemplatematchin2.c
//  rstools
//
//  Created by André Hoffmann on 4/13/14.
//
//

#include <stdio.h>
#include <string.h>
#include "rsmathutils.h"
#include "rsniftiutils.h"

void rsTemplateMatchinPrintHelp() {
    printf(
		RSTOOLS_VERSION_LABEL "\n\n"
        "Takes in a list of niftis with one ore more volumes\n"
        "which contain t-values. Those are then compared to a\n"
        "template by first thresholding all volumes to the\n"
        "specified t-value and then computing precision and\n"
        "recall based on the number of voxels exceeding/not\n"
        "meeting the threshold in both the template and the\n"
        "tested volume.\n"
        "basic usage:  stdin | rstemplatematching2 -output <volume> -n <int>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -precision <txt>       : the txt file in which a matrix with the\n"
        "                            precision values will be saved\n"
    );

	printf(
    	"   -recall <txt>          : the txt file in which a matrix with the\n"
    	"                            recall values will be saved\n"
	);
	
	printf(
    	"   -mask <nii>            : a mask specifying the region of interest\n"
    	"                            on which the computations will be\n"
        "                            performed on\n"
	);
    
	printf(
    	"   -template <nii>        : template file that will be compared against\n"
	);

	printf(
    	"   -q <float>             : FDR q-limit that will be used as a threshold\n"
	);

    printf(
        "   -threads <int>         : number of threads used for processing\n"
    );
    
    printf(
        "   -v[erbose]             : show debug information\n"
        "\n"
    );
}

struct rsInputFile {
    char*   path;
	BOOL    readable;
    FSLIO*  fslio;
    short   xDim;
    short   yDim;
    short   zDim;
    short   vDim;
    size_t  dt;
    short   pixtype;
    float   inter;
    float   slope;
    double* data;
	float   dof;
	double  tThreshold;
};


struct rsInputFile rsOpenInputFile(char* path) {
    struct rsInputFile f;
    
    f.path     = (char*)malloc((strlen(path)+1)*sizeof(char));
    strcpy(f.path, path);
    f.path     = rsTrimString(f.path);
    f.readable = FALSE;
    
    /* open file */
    f.fslio = FslOpen(f.path, "rb");
    
    if (f.fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n", f.path);
        return f;
    }
    
	/* determine degrees of freedom */
	float unused1, unused2;
	short intent_code;
	FslGetIntent(f.fslio, &intent_code, &f.dof, &unused1, &unused2);
	
    if (intent_code != NIFTI_INTENT_TTEST) {
        fprintf(stderr, "\nFile %s does not seem to contain T-values!\n", f.path);
        return f;
    }

    f.readable = TRUE;
    
	/* determine dimensions */
	FslGetDim(f.fslio, &f.xDim, &f.yDim, &f.zDim, &f.vDim);
    
    /* determine scaling */
    f.inter = 0.0;
    f.slope = 1.0;
    
    if (f.fslio->niftiptr->scl_slope != 0) {
        f.slope = f.fslio->niftiptr->scl_slope;
        f.inter = f.fslio->niftiptr->scl_inter;
    }
	
	/* determine datatype */
	f.dt = FslGetDataType(f.fslio, &f.pixtype);


    /* read in file */
    f.data = malloc((size_t)f.xDim*(size_t)f.yDim*(size_t)f.zDim*(size_t)f.vDim*(size_t)f.dt/(size_t)8);
    
    FslReadVolumes(f.fslio, f.data, f.vDim);
    
    return f;
}

void rsCloseInputFile(struct rsInputFile* f) {
    free((*f).data);
    FslClose((*f).fslio);
    free((*f).fslio);
    (*f).readable = FALSE;
}

struct rsInputFile *rsReadFileListFromStandardInput(unsigned int *nFiles) {
    char *line = NULL;
    size_t len = 0;
    size_t read;
    int sizeFilesBuffer = 1;
    struct rsInputFile *files = (struct rsInputFile*)malloc(sizeof(struct rsInputFile)*sizeFilesBuffer);
    *nFiles=0;
    
    while ((read = getline(&line, &len, stdin)) != -1) {
        *nFiles = *nFiles + 1;
        
        // Check if we're running out of memory and extend the array if necessary
        if ( *nFiles >= sizeFilesBuffer ) {
            sizeFilesBuffer = sizeFilesBuffer + 10;
            struct rsInputFile* tmpFiles = (struct rsInputFile*)realloc(files, sizeFilesBuffer * sizeof(struct rsInputFile));
            if (tmpFiles) {
                files = tmpFiles;
            } else {
                fprintf(stderr, "Could not allocate enough memory to read the file list from stdin.\n");
                exit(EXIT_FAILURE);
            }
        }
        
        files[*nFiles-1] = rsOpenInputFile(line);
    }
    
    if (line) free(line);
    
    return files;
}

int main(int argc, char * argv[]) {
    	
	char *precisionpath = NULL;
	char *recallpath    = NULL;
	char *maskpath      = NULL;
	char *templatepath  = NULL;
	
	BOOL verbose = FALSE;
    int threads = 1;
	float threshold = -1.0;
	
	int ac;
    
	if( argc < 2 ) {
        rsTemplateMatchinPrintHelp();
        return 1;
    }
    
	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			rsTemplateMatchinPrintHelp();
            return 1;
		} else if ( ! strncmp(argv[ac], "-v", 2) ) {
			verbose = TRUE;
		} else if ( ! strcmp(argv[ac], "-precision") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -precision\n");
				return 1;
			}
			precisionpath = argv[ac];
		} else if ( ! strcmp(argv[ac], "-recall") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -recall\n");
				return 1;
			}
			recallpath = argv[ac];
		} else if ( ! strcmp(argv[ac], "-mask") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -mask\n");
				return 1;
			}
			maskpath = argv[ac];
		} else if ( ! strcmp(argv[ac], "-template") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -template\n");
				return 1;
			}
			templatepath = argv[ac];
		} else if ( ! strcmp(argv[ac], "-threads") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -threads\n");
           		return 1;
           	}
           	threads = atoi(argv[ac]);
        } else if ( ! strcmp(argv[ac], "-q") ) {
	  		if( ++ac >= argc ) {
	       		fprintf(stderr, "** missing argument for -q\n");
	       		return 1;
	       	}
	       	threshold = atof(argv[ac]);
	    } else {
			fprintf(stderr, "\nError, unrecognized command %s\n", argv[ac]);
		}
	}
    
	if ( maskpath == NULL ) {
		fprintf(stderr, "No mask specified!\n");
		return 1;
	}
	
	if ( templatepath == NULL ) {
		fprintf(stderr, "No template file specified!\n");
		return 1;
	}
	
	if ( threshold <= 0 ) {
		threshold = 0.05;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Precision file: %s\n", precisionpath);
	    fprintf(stdout, "Recall file: %s\n", recallpath);
		fprintf(stdout, "Mask file: %s\n", maskpath);
		fprintf(stdout, "Template file: %s\n", templatepath);
		fprintf(stdout, "Threshold: %.2f\n", threshold);
    }

    // Load list of files
    unsigned int nFiles = 0;
    struct rsInputFile *files = rsReadFileListFromStandardInput(&nFiles);
	size_t fileListLength = 1;
	
    for (int n=0; n<nFiles; n=n+1) {
        const struct rsInputFile file = files[n];
        
        if ( ! file.readable ) {
            fprintf(stderr, "File '%s' is not accessible.\n", file.path);
            return 1;
        }
        
        if (verbose) {
            fprintf(stdout, "File: %s, Volumes: %d\n", file.path, file.vDim);
        }

		fileListLength = fileListLength + strlen(file.path) + 5;
    }
    
    if ( nFiles < 1 ) {
        fprintf(stderr, "No files were supplied via standard input!\n");
        return 1;
    }

	char fileList[fileListLength];
	size_t bytesWritten = 0;
    for (int n=0; n<nFiles; n=n+1) {
        const struct rsInputFile file = files[n];
		sprintf(&fileList[bytesWritten], "%s,%04d\n", file.path, file.vDim);
		bytesWritten = bytesWritten + strlen(file.path) + 6;
	}
	fileList[fileListLength-1] = '\0';
    
    const struct rsInputFile refFile = files[0];
	struct rsInputFile template = rsOpenInputFile(templatepath);
	
	if ( !template.readable ) {
		fprintf(stderr, "Template file could not be read!\n");
        return 1;
	}
    
	// Load mask
    unsigned long nPoints = 0L;
    Point3D *maskPoints = ReadMask(maskpath, refFile.xDim, refFile.yDim, refFile.zDim, &nPoints, NULL, refFile.fslio, NULL);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        return 1;
    }

	// Find the file with the largest DOF
	struct rsInputFile fileWithLargestDOF = files[0];
	for ( int row = 1; row < nFiles; row = row + 1 ) {
		const struct rsInputFile file = files[row];
		if ( file.dof > fileWithLargestDOF.dof ) {
			fileWithLargestDOF = file;
		}
	}
	
	// Compute FDR-threshold for the file with the largest DOF
	struct rsFDRResult fdrResult;
	double ***data;
	data = d3matrix(fileWithLargestDOF.zDim-1,  fileWithLargestDOF.yDim-1, fileWithLargestDOF.xDim-1);
	rsExtractVolumeFromBuffer(fileWithLargestDOF.fslio, data[0][0], fileWithLargestDOF.data, fileWithLargestDOF.slope, fileWithLargestDOF.inter, fileWithLargestDOF.vDim-1, fileWithLargestDOF.xDim, fileWithLargestDOF.yDim, fileWithLargestDOF.zDim);
	fdrResult = rsComputeTThresholdFDR(data, threshold, maskPoints, nPoints, fileWithLargestDOF.dof);
	
	// Compute corresponding t-threshold for all other files
	for ( int row = 0; row < nFiles; row = row + 1 ) {
		struct rsInputFile file = files[row];
		file.tThreshold = rsComputeTValueFromPValue(fdrResult.p / 2, file.dof);
		if (verbose) {
			fprintf(stdout, "File: %s, DOF: %02.0f, t-FDR: %02.9f, uncorrected-p: %02.9f\n", file.path, file.dof, file.tThreshold, fdrResult.p);
		}
	}

	// Prepare processing
//	gsl_rng *randgen = rsGetRandomNumberGenerator();
//	double *tValues;
//	double *series;
//	size_t *indices;
//	size_t fIndex;
//    short t,x,y,z,f;
//    
//	for ( int row = 0; row < nFiles; row = row + 1 ) {
//		const struct rsInputFile file = files[row];
//		
//		for ( int col = 0; col < file.vDim; col = col + 1 ) {
//			double ***data = d3matrix(file.zDim-1,  file.yDim-1, file.xDim-1);
//			rsExtractVolumeFromBuffer(file.fslio, data[0][0], file.data, file.slope, file.inter, threshold, file.xDim, file.yDim, file.zDim);
//			
//			struct rsFDRResult fdrResult = rsComputeTThresholdFDR(data, threshold, maskPoints, nPoints, file.dof);
//			
//			if (verbose) {
//				fprintf(stdout, "Row: % 2d, Col: %02d, DOF: %02.0f, t-FDR: %02.9f, p-FDR: %02.9f, i-FDR: %d, iNorm-FDR: %.4f\n", row, col, file.dof, fdrResult.t, fdrResult.p, fdrResult.i, fdrResult.iNormalized);
//			}
//			
//			free(data[0][0]);
//			free(data[0]);
//			free(data);
//		}
//	}

//    #pragma omp parallel num_threads(threads) private(x,y,t,f,tValues,indices,fIndex,series) shared(randgen,buffer,nFiles,files,nOutputVolumes, nNiftis)
//    {
//        #pragma omp for schedule(guided)
//        for (z=0; z<refFile.zDim; z=z+1) {
//			
//			tValues = malloc((size_t)nOutputVolumes * sizeof(double));
//            indices = malloc((size_t)nFiles         * sizeof(size_t));
//			series  = malloc((size_t)nNiftis        * sizeof(double));
//
//			// prepare array with indices that will be shuffled
//			for (size_t i = 0; i < (size_t)nFiles; i=i+1) {
//				indices[i] = i;
//			}
//	
//            for (y=0; y<refFile.yDim; y=y+1) {
//                for (x=0; x<refFile.xDim; x=x+1) {
//                    Point3D point = MakePoint3D(x, y, z);
//
//					// shuffle timepoints
//					gsl_ran_shuffle(randgen, &indices[0], nFiles, sizeof(size_t));
//
//                    for (t=0; t<refFile.vDim; t=t+1) {
//                        for (f=0; f<(int)nNiftis; f=f+1) {
//							fIndex = indices[f];
//                            const struct rsInputFile *file = &files[fIndex];
//
//                            rsExtractPointsFromBuffer(
//                                (*file).fslio,
//                                &series[f],
//                                (*file).data,
//                                (*file).slope,
//                                (*file).inter,
//                                &point,
//                                1L,
//                                t,
//                                (*file).xDim,
//                                (*file).yDim,
//                                (*file).zDim,
//                                (*file).vDim
//                            );
//                        }
//                        
//                        tValues[t] = rsOneSampleTTest(series, nNiftis, 0.0);
//                    }   
//                    
//                    rsWriteTimecourseToBuffer(fslioOutput, tValues, buffer, refFile.slope, refFile.inter, point, refFile.xDim, refFile.yDim, refFile.zDim, nOutputVolumes);
//                }
//            }
//			
//			free(series);
//			free(indices);
//            free(tValues);
//        }
//    }
//    
//	rsDestroyRandomNumberGenerator();
//    FslWriteVolumes(fslioOutput, buffer, nOutputVolumes);
    
    // Close files
    
    for (int n=0; n<nFiles; n=n+1) {
        struct rsInputFile file = files[n];
        rsCloseInputFile(&file);
    }
    rsCloseInputFile(&template);
    free(files);
}