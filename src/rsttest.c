#include <stdio.h>
#include <string.h>
#include "maths/rsmathutils.h"
#include "nifti/rsniftiutils.h"

void rsTTestPrintHelp() {
    printf(
        RSTOOLS_VERSION_LABEL "\n\n"
        "Takes in a list of niftis with one ore more volumes\n"
        "or a single nifti with more than one volume via stdin\n"
        "and performs a one-sample t-test on it.\n"
        "In the case of multiple niftis, volumes with the same\n"
        "index will be compared. If only a single nifti is supplied\n"
        "all of its volumes will be compared to each other.\n"
        "\n"
        "basic usage:  stdin | rsttest -output <volume>\n"
        "\n"
    );
    
    printf(
        "options:\n"
    );
    
    printf(
        "   -output <volume>       : the volume in which the result will be saved\n"
    );
    
    printf(
        "   -threads <int>         : number of threads used for processing\n"
    );
    
    printf(
        "   -v[erbose]             : show debug information\n"
        "\n"
    );
}

rsNiftiFile **rsReadFileListFromStandardInput(unsigned int *nFiles) {
    char *line = NULL;
    size_t len = 0;
    size_t read;
    int sizeFilesBuffer = 1;
    rsNiftiFile **files = (rsNiftiFile**)rsMalloc(sizeof(rsNiftiFile*)*sizeFilesBuffer);
    *nFiles=0;
    
    while ((read = getline(&line, &len, stdin)) != -1) {
        *nFiles = *nFiles + 1;
        
        // Check if we're running out of memory and extend the array if necessary
        if ( *nFiles >= sizeFilesBuffer ) {
            sizeFilesBuffer = sizeFilesBuffer + 10;
            rsNiftiFile **tmpFiles = (rsNiftiFile**)realloc(files, sizeFilesBuffer * sizeof(rsNiftiFile*));
            if (tmpFiles) {
                files = tmpFiles;
            } else {
                fprintf(stderr, "Could not allocate enough memory to read the file list from stdin.\n");
                exit(EXIT_FAILURE);
            }
        }
        
        files[*nFiles-1] = rsOpenNiftiFile(line, RSNIFTI_OPEN_READ);
    }
    
    if (line) free(line);
    
    return files;
}

int main(int argc, char * argv[]) {
    
    char *outputpath = NULL;
    
    BOOL verbose = FALSE;
    int threads = 1;
    
    int ac;
    
    if( argc < 2 ) {
        rsTTestPrintHelp();
        return 1;
    }
    
    /* parse parameters */
    for( ac = 1; ac < argc; ac++ ) {
        if( ! strncmp(argv[ac], "-h", 2) ) {
            rsTTestPrintHelp();
            return 1;
        } else if ( ! strncmp(argv[ac], "-v", 2) ) {
            verbose = TRUE;
        } else if ( ! strcmp(argv[ac], "-output") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -output\n");
                return 1;
            }
            outputpath = argv[ac];  /* no string copy, just pointer assignment */
        } else if ( ! strcmp(argv[ac], "-threads") ) {
            if( ++ac >= argc ) {
                fprintf(stderr, "** missing argument for -threads\n");
                return 1;
            }
            threads = atoi(argv[ac]);  /* no string copy, just pointer assignment */
        } else {
            fprintf(stderr, "\nError, unrecognized command %s\n", argv[ac]);
        }
    }
    
    if ( outputpath == NULL ) {
        fprintf(stderr, "No output volume specified!\n");
        return 1;
    }
    
    if ( verbose ) {
        fprintf(stdout, "Output file: %s\n", outputpath);
    }
    
    // Load list of files
    unsigned int nFiles = 0;
    rsNiftiFile **files = rsReadFileListFromStandardInput(&nFiles);
    size_t fileListLength = 1;
    
    for (int n=0; n<nFiles; n=n+1) {
        const rsNiftiFile *file = files[n];
        
        if ( ! file->readable ) {
            fprintf(stderr, "File '%s' is not accessible.\n", file->path);
            return 1;
        }
        
        if (verbose) {
            fprintf(stdout, "File: %s, Volumes: %d\n", file->path, file->vDim);
        }

        fileListLength = fileListLength + strlen(file->path) + 2 + (size_t)fmaxf(rsCountDigits(file->vDim), 4);
    }
    
    if ( nFiles < 1 ) {
        fprintf(stderr, "No files were supplied via standard input!\n");
        return 1;
    }

    // Prepare comment containing the file list
    char fileList[fileListLength];
    size_t bytesWritten = 0;
    for (int n=0; n<nFiles; n=n+1) {
        const rsNiftiFile *file = files[n];
        sprintf(&fileList[bytesWritten], "%s,%04d\n", file->path, file->vDim);
        bytesWritten = bytesWritten + strlen(file->path) + 2 + (size_t)fmaxf(rsCountDigits(file->vDim), 4);
    }
    fileList[fileListLength-1] = '\0';
    
    char *comment1 = rsMergeStringArray(argc, argv);
    char *comment2 = "\nFilelist:\n";
    size_t commentLength = strlen(comment1)+strlen(comment2)+fileListLength+1;
    char comment[commentLength];
    sprintf(&comment[0], "%s%s%s\n", comment1, comment2, fileList);
    comment[commentLength-1] = '\0';
    
    // Prepare output file
    const rsNiftiFile *refFile = files[0];
    const size_t nOutputVolumes = (nFiles > 1) ? refFile->vDim : 1;
    rsNiftiFile *outputFile = rsCloneNiftiFile(outputpath, refFile, RSNIFTI_OPEN_ALLOC, nOutputVolumes);
    
    if ( ! outputFile->readable ) {
        exit(EXIT_FAILURE);
    }
    
    FslSetIntent(outputFile->fslio, NIFTI_INTENT_TTEST, (nFiles > 1) ? (nFiles-1) : (refFile->vDim-1), 0, 0);
    rsWriteNiftiHeader(outputFile->fslio, &comment[0]);
        
    short t,x,y,z;
    Point3D *point;
    
    // Iterate over all voxels in the nifti
    #pragma omp parallel num_threads(threads) private(x,y,t,point) shared(outputFile,nFiles,files)
    {
        #pragma omp for schedule(guided)
        for (z=0; z<refFile->zDim; z++) {
            
            for (y=0; y<refFile->yDim; y=y+1) {
                for (x=0; x<refFile->xDim; x=x+1) {
                    point = rsMakePoint3D(x, y, z);
                    
                    double *tValues = (double*)malloc(sizeof(double)*nOutputVolumes);
                    
                    // if more than one file is supplied compare the same volume indices with each other
                    if ( nFiles > 1 ) {
                        for (t=0; t<refFile->vDim; t=t+1) {
                            double *series = (double*)rsMalloc(sizeof(double)*nFiles);

                            for (int f=0; f<nFiles; f=f+1) {
                                const rsNiftiFile *file = files[f];
                                rsExtractPointsFromBuffer(
                                    file->dt,
                                    &series[f],
                                    file->data,
                                    file->slope,
                                    file->inter,
                                    point,
                                    1L,
                                    t,
                                    file->xDim,
                                    file->yDim,
                                    file->zDim,
                                    file->vDim
                                );
                            }
                            
                            tValues[t] = rsOneSampleTTest(series, nFiles, 0.0);
                            free(series);
                        }                            
                    } else { // otherwise perform the t-test along the different volumes in the 4D-nifti
                        double *series = (double*)malloc(sizeof(double)*refFile->vDim);
                        
                        for (t=0; t<refFile->vDim; t=t+1) {
                            
                            rsExtractPointsFromBuffer(
                                refFile->dt,
                                &series[t],
                                refFile->data,
                                refFile->slope,
                                refFile->inter,
                                point,
                                1L,
                                t,
                                refFile->xDim,
                                refFile->yDim,
                                refFile->zDim,
                                refFile->vDim
                            );
                        }
                        
                        tValues[0] = rsOneSampleTTest(series, refFile->vDim, 0.0);
                        free(series);                            
                    }
                    
                    rsWriteTimecourseToBuffer(outputFile->dt, tValues, outputFile->data, refFile->slope, refFile->inter, point, refFile->xDim, refFile->yDim, refFile->zDim, nOutputVolumes);
                    free(tValues);
                    free(point);
                }
            }
        }
    }
    
    // Write result
    
    FslWriteVolumes(outputFile->fslio, outputFile->data, nOutputVolumes);
    
    // Close files
    
    for (int n=0; n<nFiles; n=n+1) {
        rsNiftiFile *file = files[n];
        rsCloseNiftiFileAndFree(file);
    }
    
    rsCloseNiftiFileAndFree(outputFile);
    free(files);
}
