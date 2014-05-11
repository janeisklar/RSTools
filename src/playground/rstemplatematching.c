//
//  rstemplatematchin2.c
//  rstools
//
//  Created by André Hoffmann on 4/13/14.
//
//

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "src/maths/rsmathutils.h"
#include "src/nifti/rsniftiutils.h"
#include "matio.h"

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
		"For all output matrices the files will be in the row\n"
		"and the volumes in the columns."
        "basic usage:  stdin | rstemplatematching2 -output <volume> -n <int>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -output <mat>          : a mat file in which all output variables\n"
        "                            will be stored\n"
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
    	"   -pTemplate <float>     : uncorrected p-threshold that's used for the\n"
        "                            template\n"
	);
	
	printf(
    	"   -steps <float>         : number of thresholds that should be tested\n"
        "                            within the supplied range of p-values\n"
	);

	printf(
		"   -p <float>             : range of uncorrected p-values that will be\n"
		"                            used to threshold both the input volumes\n"
		"                            and the template before comparing them\n"
		"                            (defaults to p<0.001)\n"
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
    	
	char *outputpath    = NULL;
	char *maskpath      = NULL;
	char *templatepath  = NULL;
	
	BOOL verbose = FALSE;
    int threads = 1;
	double pRange[2] = {0.0005, 0.01};
	size_t pSteps = 95;
	double pTemplate = 0.001;
	
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
		} else if ( ! strcmp(argv[ac], "-output") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -output\n");
				return 1;
			}
			outputpath = argv[ac];
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
		} else if ( ! strcmp(argv[ac], "-pTemplate") ) {
  			if( ++ac >= argc ) {
           		fprintf(stderr, "** missing argument for -pTemplate\n");
           		return 1;
           	}
           	pTemplate = atof(argv[ac]);
        } else if ( ! strcmp(argv[ac], "-steps") ) {
 			if( ++ac >= argc ) {
          		fprintf(stderr, "** missing argument for -steps\n");
          		return 1;
          	}
          	pSteps = atoi(argv[ac]);
        } else if ( ! strcmp(argv[ac], "-p") ) {
	  		if( (++ac+1) >= argc ) {
	       		fprintf(stderr, "** missing argument(s) for -p\n");
	       		return 1;
	       	}
	       	pRange[0] = atof(argv[ac]);
			++ac;
	       	pRange[1] = atof(argv[ac]);			
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
	
	if ( pRange[0] <= 0 || pRange[1] <= 0 || pRange[0] > pRange[1] ) {
		fprintf(stderr, "Given thresholds are invalid!\n");
		return 1;
	}
	
    if ( verbose ) {
		fprintf(stdout, "Mask file: %s\n", maskpath);
		fprintf(stdout, "Template file: %s\n", templatepath);
		fprintf(stdout, "p-threshold range: %.9f...%.9f\n", pRange[0], pRange[1]);
		fprintf(stdout, "Template p-threshold: %.9f\n", pTemplate);
		fprintf(stdout, "Number of p-values to test: %ld\n", pSteps);
		fprintf(stdout, "Writing statistics to: %s\n", outputpath);
    }

	rsSetThreadsNum(threads);
	
	if ( verbose ) {
		fprintf(stdout, "\nReading in nifti-files..\n");
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

		fileListLength = fileListLength + strlen(file.path) + 2 + (size_t)fmaxf(rsCountDigits(file.vDim), 4);
    }


	// Open output file
	mat_t *mat;
	matvar_t *struct_fields[10] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
	mat = Mat_CreateVer(outputpath,NULL,MAT_FT_MAT73);
	
	if ( mat == NULL ) {
		fprintf(stderr, "Mat file could not be written to '%s'\n", outputpath);
		return 1;
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
		bytesWritten = bytesWritten + strlen(file.path) + 2 + (size_t)fmaxf(rsCountDigits(file.vDim), 4);
	}
	fileList[fileListLength-1] = '\0';
    
    const struct rsInputFile refFile = files[0];
	struct rsInputFile template = rsOpenInputFile(templatepath);
	
	if ( !template.readable ) {
		fprintf(stderr, "Template file could not be read!\n");
        return 1;
	}
    
	// Load mask
	double ***mask = d3matrix(refFile.zDim, refFile.yDim, refFile.xDim);
    unsigned long nPoints = 0L;
    Point3D *maskPoints = ReadMask(maskpath, refFile.xDim, refFile.yDim, refFile.zDim, &nPoints, NULL, refFile.fslio, mask);
    if ( maskPoints == NULL) {
        fprintf(stderr, "\nError: Mask invalid.\n");
        return 1;
    }
	
	if (verbose) {
		fprintf(stdout, "\nComputing thresholds\n");
	}
	
	// compute thresholds that are going to be tested
	double pThresholds[pSteps];
	pThresholds[0] = pRange[0];
	const double pRangeWidth = pRange[1]-pRange[0];
	for ( int i=1; i<pSteps; i=i+1 ) {
		// sample thresholds in a way that lower probabilites are tested more thoroughly than higher ones
		// as they will show more changes in the rate of classification
		const double x = (double)(i) / (double)(pSteps-1.0); // walk from 1/nSteps to 1
		pThresholds[i] = pRange[0] + pow(x,10) * pRangeWidth;
	}
	
	// Compute corresponding t-threshold for the template 
	template.tThreshold = rsComputeTValueFromPValue(pTemplate / 2.0, template.dof);
	if (verbose) {
		fprintf(stdout, "File: %s, DOF: %02.0f, t: %02.9f, uncorrected-p: %02.9f\n", template.path, template.dof, template.tThreshold, pTemplate);
	}
	
	// Prepare actual computation
	short x, y, z, row, col, processedThresholds = 0, nColumns=refFile.vDim;
	int p;
	double ***hits              = d3matrix(pSteps-1, nFiles-1, refFile.vDim-1); // aka true positives
	double ***misses            = d3matrix(pSteps-1, nFiles-1, refFile.vDim-1); // aka false negatives
	double ***falseAlarms       = d3matrix(pSteps-1, nFiles-1, refFile.vDim-1); // aka false positives
	double ***correctRejections = d3matrix(pSteps-1, nFiles-1, refFile.vDim-1); // aka true negatives
	double ***templateData      = d3matrix(refFile.zDim-1, refFile.yDim-1, refFile.xDim-1);
    rsExtractVolumeFromBuffer(template.fslio, templateData[0][0], template.data, template.slope, template.inter, 0, refFile.xDim, refFile.yDim, refFile.zDim);
	
	// Reset statistics
	for ( p=0; p<pSteps; p=p+1 ) {
		for (row=0; row<nFiles; row=row+1) {
			for (col=0; col<nColumns; col=col+1) {
				hits[p][row][col]              = 0;
				misses[p][row][col]            = 0;
				falseAlarms[p][row][col]       = 0;
				correctRejections[p][row][col] = 0;
			}
		}
	}
	
	if (verbose) {
		fprintf(stdout, "\nComparing volumes to the template\n");
	}
	
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(p,z,y,x,row,col) shared(processedThresholds,templateData,template,hits,misses,falseAlarms,nFiles,files,mask)
    {
		/* Iterate over all thresholds, voxels, files(rows) and volumes(columns) in the mask that should be compared to the template */
        #pragma omp for schedule(guided)
		for ( p = 0; p < pSteps; p = p + 1 ) {
 	       for ( row = 0; row < nFiles; row = row + 1 ) {
				const struct rsInputFile *file = &files[row];
		
				// Compute t-threshold that corresponds to the given p-value for all files
				double tThreshold = rsComputeTValueFromPValue(pThresholds[p] / 2, (*file).dof);

		        for (z=0; z<refFile.zDim; z=z+1) {
		            for (y=0; y<refFile.yDim; y=y+1) {
		                for (x=0; x<refFile.xDim; x=x+1) {
	
		                    /* If it's not in the mask skip it to improve the performance */
		                    if (mask != NULL && mask[z][y][x] < 0.1) {
		                        continue;
		                    }

							/* Don't compare voxels correspond to a NAN-value in the template */
							if ( templateData[z][y][x] != templateData[z][y][x] ) {
								continue;
							}
					
							const BOOL voxelShouldBeSignificant = templateData[z][y][x] >= template.tThreshold;

		                    /* read out columns for the current row */
		                    double *columns = (double*)malloc(refFile.vDim*sizeof(double));
		                    rsExtractTimecourseFromBuffer(
								(*file).fslio, 
								&columns[0], 
								(*file).data, 
								(*file).slope, 
								(*file).inter, 
								MakePoint3D(x,y,z), 
								(*file).xDim, 
								(*file).yDim, 
								(*file).zDim, 
								(*file).vDim
							);
                	
							for ( col = 0; col<refFile.vDim; col=col+1 ) {
								const BOOL voxelSignificant = columns[col] >= tThreshold;

								if ( voxelShouldBeSignificant ) {
									if ( voxelSignificant ) {
										hits[p][row][col] += 1;
									} else {
										misses[p][row][col] += 1;
									}
								} else {
									if ( voxelSignificant ) {
										falseAlarms[p][row][col] += 1;
									} else {
										correctRejections[p][row][col] += 1;
									}
								}
							}
                    
		                    free(columns);
		                }
		            }
		        }
			}
			
			/* show progress */
			if (verbose) {
            	#pragma omp atomic
            	processedThresholds += 1;
       
            	if (pSteps > 9 && processedThresholds > 0 && processedThresholds % (short)(pSteps / 10) == 0) {
                	fprintf(stdout, "..%.0f%%\n", ceil((float)processedThresholds*100.0 / (float)pSteps));
            	}
			}
		}

		// Create output mat-file
	    matvar_t *struct_matvar;
		//matvar_t *struct_fields[1] = {NULL};
	    size_t    sDims[3] = {nFiles,nColumns,pSteps};
	    int       iDims[3] = {nFiles,nColumns,pSteps};
		size_t    length = (size_t)pSteps*(size_t)nFiles*(size_t)nColumns;
		
		#pragma omp sections
		{			
			// Close files and cleanup
			#pragma omp section
			{
				for (int n=0; n<nFiles; n=n+1) {
			        struct rsInputFile file = files[n];
			        rsCloseInputFile(&file);
			    }
			
				rsCloseInputFile(&template);
			    free(files);
				free(mask);
				free(templateData);
			}
			
			// Prepare hits, misses, correctRejections, falseAlarms
			#pragma omp section
			{
				double *_hits              = malloc(length * sizeof(double));
				double *_misses            = malloc(length * sizeof(double));
				double *_correctRejections = malloc(length * sizeof(double));
				double *_falseAlarms       = malloc(length * sizeof(double));
				int *index = malloc(3 * sizeof(int));
				
				for ( p = 0; p < pSteps; p = p + 1 ) {              index[2] = p+1;
	            	for ( row = 0; row < nFiles; row = row + 1 ) {  index[0] = row+1;
						for ( col = 0; col<nColumns; col=col+1 ) {  index[1] = col+1;
							const int lindex = Mat_CalcSingleSubscript(3, iDims, index);
							_hits[lindex]              = hits[p][row][col];
							_falseAlarms[lindex]       = falseAlarms[p][row][col];
							_correctRejections[lindex] = correctRejections[p][row][col];
							_misses[lindex]            = misses[p][row][col];
						}
					}
				}
				
				struct_fields[0] = Mat_VarCreate("hits",              MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims,   _hits,              0);
				struct_fields[1] = Mat_VarCreate("misses",            MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims,   _misses,            0);
				struct_fields[2] = Mat_VarCreate("falseAlarms",       MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims,   _falseAlarms,       0);
				struct_fields[3] = Mat_VarCreate("correctRejections", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims,   _correctRejections, 0);
				struct_fields[4] = Mat_VarCreate("thresholds",        MAT_C_DOUBLE, MAT_T_DOUBLE, 1, &pSteps, &pThresholds[0],    0);

				free(index);
				free(_hits);
				free(_misses);
				free(_correctRejections);
				free(_falseAlarms);
			}
			
			// Compute Precision
			#pragma omp section
			{
				double *precision = malloc(length * sizeof(double));
				int *index = malloc(3 * sizeof(int));
				
				for ( p = 0; p < pSteps; p = p + 1 ) {              index[2] = p+1;
	            	for ( row = 0; row < nFiles; row = row + 1 ) {  index[0] = row+1;
						for ( col = 0; col<nColumns; col=col+1 ) {  index[1] = col+1;
							const int lindex = Mat_CalcSingleSubscript(3, iDims, index);
							precision[lindex] = hits[p][row][col] / (hits[p][row][col] + falseAlarms[p][row][col]);
						}
					}
				}
				struct_fields[5] = Mat_VarCreate("precision", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims, precision, 0);
				
				free(index);
				free(precision);
			}
			
			// Compute Recall
			#pragma omp section
			{
				double *recall = malloc(length * sizeof(double));
				int *index = malloc(3 * sizeof(int));
				
				for ( p = 0; p < pSteps; p = p + 1 ) {              index[2] = p+1;
	            	for ( row = 0; row < nFiles; row = row + 1 ) {  index[0] = row+1;
						for ( col = 0; col<nColumns; col=col+1 ) {  index[1] = col+1;
							const int lindex = Mat_CalcSingleSubscript(3, iDims, index);
							recall[lindex] = hits[p][row][col] / (hits[p][row][col] + misses[p][row][col]);
						}
					}
				}
				struct_fields[6] = Mat_VarCreate("recall", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims, recall, 0);

				free(index);
				free(recall);
			}
			
			// Compute F1 score
			// @see http://en.wikipedia.org/wiki/F1_score
			#pragma omp section
			{
				double *f1score = malloc(length * sizeof(double));
				int *index = malloc(3 * sizeof(int));
				
				for ( p = 0; p < pSteps; p = p + 1 ) {              index[2] = p+1;
	            	for ( row = 0; row < nFiles; row = row + 1 ) {  index[0] = row+1;
						for ( col = 0; col<nColumns; col=col+1 ) {  index[1] = col+1;
							const int lindex = Mat_CalcSingleSubscript(3, iDims, index);
							const double precision = hits[p][row][col] / (hits[p][row][col] + falseAlarms[p][row][col]);
							const double recall    = hits[p][row][col] / (hits[p][row][col] + misses[p][row][col]);
							f1score[lindex] = (precision * recall) / (precision + recall);
						}
					}
				}
				struct_fields[7] = Mat_VarCreate("f1score", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims, f1score, 0);

				free(index);
				free(f1score);
			}
			
			// Compute Matthews correlation coefficient
			// @see http://en.wikipedia.org/wiki/Matthews_correlation_coefficient
			#pragma omp section
			{
				double *mcc = malloc(length * sizeof(double));
				int *index = malloc(3 * sizeof(int));
				
				for ( p = 0; p < pSteps; p = p + 1 ) {              index[2] = p+1;
	            	for ( row = 0; row < nFiles; row = row + 1 ) {  index[0] = row+1;
						for ( col = 0; col<nColumns; col=col+1 ) {  index[1] = col+1;
							const int lindex = Mat_CalcSingleSubscript(3, iDims, index);
							const double enumerator  = hits[p][row][col]*correctRejections[p][row][col] - falseAlarms[p][row][col]*misses[p][row][col];
							const double denominator = sqrt((hits[p][row][col]+falseAlarms[p][row][col])*(hits[p][row][col]+misses[p][row][col])*(correctRejections[p][row][col]+falseAlarms[p][row][col])*(correctRejections[p][row][col]+misses[p][row][col]));
							mcc[lindex] = enumerator / denominator;
						}
					}
				}
				struct_fields[8] = Mat_VarCreate("mcc", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims, mcc, 0);

				free(index);
				free(mcc);
			}
			
			// Compute Accuracy
			// @see http://en.wikipedia.org/wiki/Accuracy
			#pragma omp section
			{
				double *acc = malloc(length * sizeof(double));
				int *index = malloc(3 * sizeof(int));
				
				for ( p = 0; p < pSteps; p = p + 1 ) {              index[2] = p+1;
	            	for ( row = 0; row < nFiles; row = row + 1 ) {  index[0] = row+1;
						for ( col = 0; col<nColumns; col=col+1 ) {  index[1] = col+1;
							const int lindex = Mat_CalcSingleSubscript(3, iDims, index);
							const double enumerator  = hits[p][row][col] + correctRejections[p][row][col];
							const double denominator = enumerator +falseAlarms[p][row][col] + misses[p][row][col];
							acc[lindex] = enumerator / denominator;
						}
					}
				}
				struct_fields[9] = Mat_VarCreate("accuracy", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, sDims, acc, 0);

				free(index);
				free(acc);
			}
		}		
    }

	for ( int i=0; i<10; i=i+1 ) {
		Mat_VarWrite(mat, struct_fields[i], MAT_COMPRESSION_NONE);
		Mat_VarFree(struct_fields[i]);
	}

	Mat_Close(mat);
	free(hits);
	free(misses);
	free(falseAlarms);
	free(correctRejections);
}
