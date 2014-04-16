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
		"For all output matrices the files will be in the row\n"
		"and the volumes in the columns."
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
    	"   -f1score <txt>         : the txt file in which a matrix with the\n"
    	"                            f1-scores will be saved\n"
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
		"   -p <float>             : uncorrected p-value that will be used to\n"
		"                            threshold both the input volumes and the\n"
		"                            template before comparing them\n"
		"                            (defaults to p<0.001)\n"
	);

	printf(
    	"   -q <float>             : FDR q-limit that will be used as a threshold\n"
        "                            instead of a p-value\n"
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
	char *f1scorepath   = NULL;
	char *maskpath      = NULL;
	char *templatepath  = NULL;
	
	BOOL verbose = FALSE;
    int threads = 1;
	float fdrThreshold = -1.0;
	float pThreshold = 0.001;
	
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
		} else if ( ! strcmp(argv[ac], "-f1score") ) {
            if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -f1score\n");
				return 1;
			}
			f1scorepath = argv[ac];
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
	       	fdrThreshold = atof(argv[ac]);
	    } else if ( ! strcmp(argv[ac], "-p") ) {
	  		if( ++ac >= argc ) {
	       		fprintf(stderr, "** missing argument for -p\n");
	       		return 1;
	       	}
	       	pThreshold = atof(argv[ac]);
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
	
	if ( fdrThreshold <= 0 && pThreshold <= 0 ) {
		fprintf(stderr, "Given thresholds are invalid!\n");
		return 1;
	}
	
    if ( verbose ) {
        fprintf(stdout, "Precision file: %s\n", precisionpath);
	    fprintf(stdout, "Recall file: %s\n", recallpath);
		fprintf(stdout, "Mask file: %s\n", maskpath);
		fprintf(stdout, "Template file: %s\n", templatepath);
		if ( fdrThreshold > 0 ) {
			fprintf(stdout, "FDR-threshold: %.9f\n", fdrThreshold);
		} else {
			fprintf(stdout, "p-threshold: %.9f\n", pThreshold);
		}
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
	
	// Compute FDR-threshold if requested
	if ( fdrThreshold > 0 ) {
	
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
		fdrResult = rsComputeTThresholdFDR(data, fdrThreshold, maskPoints, nPoints, fileWithLargestDOF.dof);
	
		// Compute corresponding t-threshold for all other files
		for ( int row = 0; row < nFiles; row = row + 1 ) {
			struct rsInputFile file = files[row];
			file.tThreshold = rsComputeTValueFromPValue(fdrResult.p / 2, file.dof);
			if (verbose) {
				fprintf(stdout, "File: %s, DOF: %02.0f, t-FDR: %02.9f, uncorrected-p: %02.9f\n", file.path, file.dof, file.tThreshold, fdrResult.p);
			}
		}
		
		// Compute corresponding t-threshold for the template itself
		template.tThreshold = rsComputeTValueFromPValue(fdrResult.p / 2, template.dof);
		if (verbose) {
			fprintf(stdout, "File: %s, DOF: %02.0f, t-FDR: %02.9f, uncorrected-p: %02.9f\n", template.path, template.dof, template.tThreshold, fdrResult.p);
		}
		
	} else { // use the supplied/default p-value instead
		
		// Compute t-threshold that corresponds to the given p-value for all files
		for ( int row = 0; row < nFiles; row = row + 1 ) {
			struct rsInputFile file = files[row];
			file.tThreshold = rsComputeTValueFromPValue(pThreshold / 2, file.dof);
			if (verbose) {
				fprintf(stdout, "File: %s, DOF: %02.0f, t-threshold: %02.9f, uncorrected-p: %02.9f\n", file.path, file.dof, file.tThreshold, pThreshold);
			}
		}
		
		// Compute corresponding t-threshold for the template itself
		template.tThreshold = rsComputeTValueFromPValue(pThreshold / 2, template.dof);
		if (verbose) {
			fprintf(stdout, "File: %s, DOF: %02.0f, t-FDR: %02.9f, uncorrected-p: %02.9f\n", template.path, template.dof, template.tThreshold, pThreshold);
		}
	}

	// Prepare actual computation
	short x, y, z, row, col, threadId, processedSlices = 0, nColumns=refFile.vDim;
	double **hits          = d2matrix(nFiles-1, refFile.vDim-1); // aka true positives
	double **misses        = d2matrix(nFiles-1, refFile.vDim-1); // aka false negatives
	double **falseAlarms   = d2matrix(nFiles-1, refFile.vDim-1); // aka false positives
	double ***_hits        = d3matrix(rsGetThreadsNum()-1, nFiles-1, refFile.vDim-1); // copy for multithreading
	double ***_misses      = d3matrix(rsGetThreadsNum()-1, nFiles-1, refFile.vDim-1);
	double ***_falseAlarms = d3matrix(rsGetThreadsNum()-1, nFiles-1, refFile.vDim-1);
	double ***templateData = d3matrix(refFile.zDim-1, refFile.yDim-1, refFile.xDim-1);
    rsExtractVolumeFromBuffer(template.fslio, templateData[0][0], template.data, template.slope, template.inter, 0, refFile.xDim, refFile.yDim, refFile.zDim);
	
	// Reset statistics
	for (row=0; row<nFiles; row=row+1) {
		for (col=0; col<nColumns; col=col+1) {
			
			for ( int thread=0; thread<rsGetThreadsNum(); thread=thread+1 ) {
				_hits[thread][row][col]        = 0;
				_misses[thread][row][col]      = 0;
				_falseAlarms[thread][row][col] = 0;
			}
			
			hits[row][col]        = 0;
			misses[row][col]      = 0;
			falseAlarms[row][col] = 0;
		}
	}
	
	if (verbose) {
		fprintf(stdout, "\nComparing volumes to the template\n");
	}
	
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(y,x,row,col,threadId) shared(processedSlices,templateData,template,hits,misses,falseAlarms,_hits,_misses,_falseAlarms,nFiles,files,mask)
    {
		/* Iterate over all voxels in the mask that should be compared */
        #pragma omp for schedule(guided)
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
					
					threadId = omp_get_thread_num();
					
					const BOOL voxelShouldBeSignificant = templateData[z][y][x] >= template.tThreshold;

					/* Otherwise go through the list of files and compare them to the template */
                    for ( row = 0; row < nFiles; row = row + 1 ) {
						const struct rsInputFile *file = &files[row];
						
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
							const BOOL voxelSignificant = columns[col] >= (*file).tThreshold;

							if ( voxelShouldBeSignificant ) {
								if ( voxelSignificant ) {
									_hits[threadId][row][col] += 1;
								} else {
									_misses[threadId][row][col] += 1;
								}
							} else {
								if ( voxelSignificant ) {
									_falseAlarms[threadId][row][col] += 1;
								}
							}
						}
	                    
	                    free(columns);
					}
                }
            }
			
			/* show progress */
			if (verbose) {
            	#pragma omp atomic
            	processedSlices += 1;
            
            	if (processedSlices > 0 && processedSlices % (short)(refFile.zDim / 10) == 0) {
                	fprintf(stdout, "..%.0f%%\n", ceil((float)processedSlices*100.0 / (float)refFile.zDim));
            	}
			}
        }

		
		#pragma omp sections
		{
			// Aggregate hits from all threads
			#pragma omp section
			{
				for ( int thread=0; thread<rsGetThreadsNum(); thread=thread+1 ) {
					for (row=0; row<nFiles; row=row+1) {
						for (col=0; col<nColumns; col=col+1) {
							hits[row][col] += _hits[thread][row][col];
						}
					}
				}
			}
			
			// Aggregate misses from all threads
			#pragma omp section
			{
				for ( int thread=0; thread<rsGetThreadsNum(); thread=thread+1 ) {
					for (row=0; row<nFiles; row=row+1) {
						for (col=0; col<nColumns; col=col+1) {
							misses[row][col] += _misses[thread][row][col];
						}
					}
				}
			}
			
			// Aggregate falseAlarm from all threads
			#pragma omp section
			{
				for ( int thread=0; thread<rsGetThreadsNum(); thread=thread+1 ) {
					for (row=0; row<nFiles; row=row+1) {
						for (col=0; col<nColumns; col=col+1) {
							falseAlarms[row][col] += _falseAlarms[thread][row][col];
						}
					}
				}
			}
			
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
		}

		#pragma omp sections
		{		
			// Compute Precision
			#pragma omp section
			{
				if ( precisionpath != NULL ) {

					FILE *hFile; 
					hFile = fopen(precisionpath,"w");
				
	                for ( row = 0; row < nFiles; row = row + 1 ) {
						for ( col = 0; col<nColumns; col=col+1 ) {
							const double precision = hits[row][col] / (hits[row][col] + falseAlarms[row][col]);
							fprintf(hFile, "%.10f\t", precision);
						}
						fprintf(hFile, "\n");
					}
					
					fclose(hFile);
				}
				
			}
			
			// Compute Recall
			#pragma omp section
			{
				if ( recallpath != NULL ) {

					FILE *hFile; 
					hFile = fopen(recallpath,"w");
				
					for ( row = 0; row < nFiles; row = row + 1 ) {
						for ( col = 0; col<nColumns; col=col+1 ) {
							const double recall = hits[row][col] / (hits[row][col] + misses[row][col]);
							fprintf(hFile, "%.10f\t", recall);
						}
						fprintf(hFile, "\n");
					}
					
					fclose(hFile);
				}

			}
			
			// Compute F1 score
			// @see http://en.wikipedia.org/wiki/F1_score
			#pragma omp section
			{
				if ( f1scorepath != NULL ) {

					FILE *hFile; 
					hFile = fopen(f1scorepath,"w");
				
					for ( row = 0; row < nFiles; row = row + 1 ) {
						for ( col = 0; col<nColumns; col=col+1 ) {
							const double precision = hits[row][col] / (hits[row][col] + falseAlarms[row][col]);
							const double recall    = hits[row][col] / (hits[row][col] + misses[row][col]);
							const double f1score   = (precision * recall) / (precision + recall);
							fprintf(hFile, "%.10f\t", f1score);
						}
						fprintf(hFile, "\n");
					}
					
					fclose(hFile);
				}
			}
		}
    }
	
	free(hits);
	free(misses);
	free(falseAlarms);
}