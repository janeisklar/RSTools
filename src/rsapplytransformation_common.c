#include <src/nifti/headerinfo.h>
#include <src/nifti/rsniftiutils.h>
#include <nifti1_io.h>
#include <externals/fslio/fslio.h>
#include "nifti/rsniftiutils.h"
#include "rsapplytransformation_common.h"
#include "utils/rsio.h"
#include "rsapplytransformation_ui.h"
#include "rszeropadding_ui.h"
#include "rszeropadding_common.h"
#include "math.h"
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_vector_double.h>

typedef struct {
    char *tmpDirPath;
    rsNiftiFile *source;
    rsNiftiFile *target;
    rsApplyTransformationTransSpecification** specs;
    size_t nTransformations;
    char *referencepath;
    short volumeIndex;
    BOOL verbose;
    char *inputSuffix;
    char *inputMaskSuffix;
    char *paddedInputSuffix;
    char *outputSuffix;
    char *outputMaskSuffix;
    char *paddedOutputSuffix;
    BOOL resourceTransformation;
    BOOL keepFiles;
    char *antsPath;
} rsApplyTransformationApplyParams;

void rsApplyTransformationLoadTransformationFile(char ***output, size_t *nTransformations, FILE *file);
BOOL rsApplyTransformationParseTransformationFile(rsApplyTransformationTransSpecification*** transformations, char **transformationFile, const size_t nTransformations);
BOOL rsApplyTransformationConvertMcFlirtTransformations(rsNiftiFile* input, char ***transformations, size_t nVolumes, const char* transformationPath);
void rsApplyTransformationConvertMcFlirtTransformMatrixToAntsTransformMatrix(mat44 *antsTransform, const rsNiftiFile* input, const mat44 *M);
BOOL rsApplyTransformationApplyToVolume(const rsApplyTransformationApplyParams *params);
char *rsApplyTransformationGetMcFlirtTransformationPath(const char *tmpDirPath, const int volumeIndex, const int transformationId);
char *rsApplyTransformationGetANTsPath();
BOOL rsApplyTransformationRunANTs(const rsApplyTransformationApplyParams *params, const char *input, const char *output, BOOL highQuality);
BOOL rsApplyTransformationConvertFugueShiftToANTsWarp(const rsNiftiFile* input, const rsNiftiFile* shift, const char* warpPath, BOOL verbose);

void rsApplyTransformationInit(rsApplyTransformationParameters *p)
{
    p->parametersValid = FALSE;

    p->antsPath = rsApplyTransformationGetANTsPath();

    if (p->antsPath == NULL) {
        return;
    }

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        (const char*)p->transformationpath,
        (const char*)p->referencepath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->outputpath,
        RSIO_LASTFILE
    });
    
    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }

    rsSetThreadsNum(p->threads);
        
    /* open input file */
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }

    p->transform = fopen(p->transformationpath, "r");

    if ( p->transform == NULL ) {
        fprintf(stderr, "\nError: The transformation file that was supplied as an input (%s) could not be read.\n", p->transformationpath);
        return;
    }

    /* output the most important parameters to the user */
    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
        fprintf(stdout, "Input Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
    }
	
    p->parametersValid = TRUE;
}

void rsApplyTransformationRun(rsApplyTransformationParameters *p)
{
    p->parametersValid = FALSE;

    // load transformations
    char **transformations = NULL;
    rsApplyTransformationLoadTransformationFile(&transformations, &p->nTransformations, p->transform);

    // parse transformations
    if (!rsApplyTransformationParseTransformationFile(&p->specs, transformations, p->nTransformations)) {
        fprintf(stderr, "\nError: The transformation file that was supplied as an input (%s) could not be parsed!\n", p->transformationpath);
        return;
    }

    if (p->verbose) {
        fprintf(stdout, "Transformations: \n");
        for (short t=0; t<p->nTransformations; t++) {
            const rsApplyTransformationTransSpecification *spec = p->specs[t];
            if (spec->type == TRANS_MCFLIRT) {
                fprintf(stdout, "[MCFLIRT] ");
            } else if (spec->type == TRANS_ANTS) {
                fprintf(stdout, "[ANTS] ");
            } else if (spec->type == TRANS_MULTIPLICATION) {
                fprintf(stdout, "[MULT] ");
            } else if (spec->type == TRANS_DIVISION) {
                fprintf(stdout, "[DIV] ");
            } else if (spec->type == TRANS_FUGUE) {
                fprintf(stdout, "[FUGUE] ");
            }
            fprintf(stdout, "%s\n", spec->file);
        }
        fprintf(stdout, "\n");
    }

    // convert transformations if necessary
    for (short t=0; t<p->nTransformations; t++) {
        if (p->specs[t]->type == TRANS_MCFLIRT) {
            if (!rsApplyTransformationConvertMcFlirtTransformations(p->input, &p->specs[t]->mcflirtAntsFiles, p->input->vDim, p->specs[t]->file)) {
                return;
            }
        }
    }

    // create ouput volume
    rsNiftiFile *reference = rsOpenNiftiFile(p->referencepath, RSNIFTI_OPEN_NONE);
    p->output = rsCloneNiftiFile(p->outputpath, reference, RSNIFTI_OPEN_ALLOC, p->input->vDim);

    // create temporary directory
    char *tmpDirPath = rsString("/tmp/rsapplytransformation.XXXXXX");

    if (mkdtemp(tmpDirPath) == NULL) {
        fprintf(stderr, "Could not create temporary directory.\n");
        return;
    }

    if (p->verbose) {
        fprintf(stdout, "Using temporary directory %s\n", tmpDirPath);
    }

    // prepare resource files (such as the bias field correction) that will need to be transformed
    // into the target world space before they can be applied
    for (short t=p->nTransformations-1; t>=0; t--) {
        rsApplyTransformationTransSpecification *spec = p->specs[t];
        if (spec->type == TRANS_MULTIPLICATION || spec->type == TRANS_DIVISION) {
            spec->transIn = rsOpenNiftiFile(spec->file, RSNIFTI_OPEN_READ);
            char *transOutPath = rsMalloc(sizeof(char) * (strlen(tmpDirPath) + 60));
            if (spec->type == TRANS_DIVISION) {
                sprintf(transOutPath, "%s/xxx_%d_div_out.nii", tmpDirPath, (int) t);
            } else {
                sprintf(transOutPath, "%s/xxx_%d_mult_out.nii", tmpDirPath, (int) t);
            }
            spec->transOut = rsCloneNiftiFile(transOutPath, reference, RSNIFTI_OPEN_ALLOC, p->output->vDim);
            rsFree(transOutPath);
        } else if (spec->type == TRANS_MCFLIRT) {
            for (short i = 0; i < p->output->vDim; i++) {
                char *transPath = rsApplyTransformationGetMcFlirtTransformationPath(tmpDirPath, i, spec->transformationId);
                FILE *trans = fopen(transPath, "wb");
                fprintf(trans, "%s", spec->mcflirtAntsFiles[i]);
                fclose(trans);
                rsFree(transPath);
            }
        } else if (spec->type == TRANS_FUGUE) {
            spec->transIn = rsOpenNiftiFile(spec->file, RSNIFTI_OPEN_READ);
            spec->transOutPath = rsMalloc(sizeof(char) * (strlen(tmpDirPath) + 60));
            sprintf(spec->transOutPath, "%s/xxx_%d_fugue_out.nii", tmpDirPath, (int) t);
            if (!rsApplyTransformationConvertFugueShiftToANTsWarp(p->input, spec->transIn, spec->transOutPath, p->verbose)) {
                return;
            }
            rsCloseNiftiFileAndFree(spec->transIn);
        }
    }
    rsCloseNiftiFileAndFree(reference);

    // transform resource files into the target world space so that they can be applied to the files
    // that are actually to be transformed at the very end
    rsApplyTransformationApplyParams *params;
    short i, t, processedVolumes = 0;
    rsApplyTransformationTransSpecification *spec;
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,t,params,spec) shared(processedVolumes)
    {
        #pragma omp for schedule(guided)
        for (i=0; i<p->input->vDim; i++) {
            params = (rsApplyTransformationApplyParams*)rsMalloc(sizeof(rsApplyTransformationApplyParams));
            params->tmpDirPath = tmpDirPath;
            params->referencepath = p->referencepath;
            params->volumeIndex = i;
            params->verbose = p->verbose && p->threads < 2;
            params->resourceTransformation = TRUE;
            params->keepFiles = p->keepFiles;
            params->antsPath = p->antsPath;

            for (t=0; t<p->nTransformations; t++) {
                spec = p->specs[t];
                if (spec->type == TRANS_MULTIPLICATION || spec->type == TRANS_DIVISION) {
                    params->source = spec->transIn;
                    params->target = spec->transOut;
                    params->nTransformations = p->nTransformations - t;
                    params->specs = &p->specs[t];
                    params->inputSuffix = rsMalloc(sizeof(char) * 15);
                    params->outputSuffix = rsMalloc(sizeof(char) * 15);
                    params->inputMaskSuffix = rsMalloc(sizeof(char) * 15);
                    params->outputMaskSuffix = rsMalloc(sizeof(char) * 15);
                    params->paddedInputSuffix= rsMalloc(sizeof(char) * 25);
                    params->paddedOutputSuffix = rsMalloc(sizeof(char) * 25);
                    const char *prefix = spec->type == TRANS_DIVISION ? "div" : "mult";
                    sprintf(params->inputSuffix, "%d_%s_in", (int)t, prefix);
                    sprintf(params->outputSuffix, "%d_%s_out", (int)t, prefix);
                    sprintf(params->inputMaskSuffix, "%d_%s_mask_in", (int)t, prefix);
                    sprintf(params->outputMaskSuffix, "%d_%s_mask_out", (int)t, prefix);
                    sprintf(params->paddedInputSuffix, "%d_%s_in_padded", (int)t, prefix);
                    sprintf(params->paddedOutputSuffix, "%d_%s_out_padded", (int)t, prefix);
                    rsApplyTransformationApplyToVolume(params);
                }
            }
            rsFree(params);

            if (p->verbose) {
                #pragma omp atomic
                processedVolumes += 1;

                if (processedVolumes > 0 && processedVolumes % (short)((p->input->vDim*2) >= 10 ? (p->input->vDim*2) / 10 : (p->input->vDim*2)) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedVolumes*100.0 / ((float)p->input->vDim*2.0f)));
                }
            }
        }
    }

    // close resource files that aren't needed anymore
    for (short t=0; t<p->nTransformations; t++) {
        rsApplyTransformationTransSpecification *spec = p->specs[t];
        if (spec->type == TRANS_MULTIPLICATION || spec->type == TRANS_DIVISION) {
            rsCloseNiftiFileAndFree(spec->transIn);
        }
    }

    // transform actual input
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,params) shared(processedVolumes)
    {
        #pragma omp for schedule(guided)
        for (i=0; i<p->input->vDim; i++) {
            params = (rsApplyTransformationApplyParams*)rsMalloc(sizeof(rsApplyTransformationApplyParams));
            params->tmpDirPath = tmpDirPath;
            params->source = p->input;
            params->target = p->output;
            params->specs = p->specs;
            params->nTransformations = p->nTransformations;
            params->referencepath = p->referencepath;
            params->volumeIndex = i;
            params->verbose = p->verbose && p->threads < 2;
            params->inputSuffix = "input";
            params->outputSuffix = "output";
            params->inputMaskSuffix = "input_mask";
            params->outputMaskSuffix = "output_mask";
            params->paddedInputSuffix = "input_padded";
            params->paddedOutputSuffix = "output_padded";
            params->resourceTransformation = FALSE;
            params->keepFiles = p->keepFiles;
            params->antsPath = p->antsPath;
            rsApplyTransformationApplyToVolume(params);
            rsFree(params);

            if (p->verbose) {
                #pragma omp atomic
                processedVolumes += 1;

                if (processedVolumes > 0 && processedVolumes % (short)((p->input->vDim*2) >= 10 ? (p->input->vDim*2) / 10 : (p->input->vDim*2)) == 0) {
                    fprintf(stdout, "..%.0f%%\n", ceil((float)processedVolumes*100.0 / ((float)p->input->vDim*2.0f)));
                }
            }
        }
    }

    // close resource files that aren't needed anymore
    for (short t=0; t<p->nTransformations; t++) {
        rsApplyTransformationTransSpecification *spec = p->specs[t];
        if (spec->type == TRANS_MULTIPLICATION || spec->type == TRANS_DIVISION) {
            if (p->keepFiles) {
                rsWriteNiftiHeader(spec->transOut->fslio, p->callString);
                FslWriteVolumes(spec->transOut->fslio, spec->transOut->data, spec->transOut->vDim);
            }
            char *tmpPath = rsString(spec->transOut->path);
            rsCloseNiftiFileAndFree(spec->transOut);

            if (!p->keepFiles) {
                unlink(tmpPath);
            }
            rsFree(tmpPath);
        } else if (spec->type == TRANS_FUGUE) {
            if (!p->keepFiles) {
                unlink(spec->transOutPath);
            }
            rsFree(spec->transOutPath);
        }
    }

    // cleanup
    if (!p->keepFiles) {
        // cleanup transformations
        for (short t=0; t<p->nTransformations; t++) {
            rsApplyTransformationTransSpecification *spec = p->specs[t];

            if (spec->type != TRANS_MCFLIRT) {
                continue;
            }

            for (short i = 0; i < p->output->vDim; i++) {
                char *transPath = rsApplyTransformationGetMcFlirtTransformationPath(tmpDirPath, i, spec->transformationId);
                unlink(transPath);
                rsFree(transPath);
            }
        }

        // remove temporary directory
        rmdir(tmpDirPath);
    }
    rsFree(tmpDirPath);

    // write output
    rsWriteNiftiHeader(p->output->fslio, p->callString);
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

    p->parametersValid = TRUE;
}

BOOL rsApplyTransformationRunANTs(const rsApplyTransformationApplyParams *params, const char *input, const char *output, BOOL highQuality)
{
    char *callString = rsMalloc(sizeof(char)*30000);
    sprintf(
        callString,
        "%santsApplyTransforms -e 3 -d 3 -i %s -o %s -r %s",
        params->antsPath,
        input,
        output,
        params->referencepath
    );

    if (highQuality) {
        rsStringAppend(callString, " --interpolation LanczosWindowedSinc");
    }

    // append transformations to the call string
    // NOTE: iterate in the reverse order as ANTs expects the transformations to be specified in that manner
    for (short i=params->nTransformations-1; i>=0; i--) {
        const rsApplyTransformationTransSpecification *spec = params->specs[i];

        if (spec->isAppliedInTargetSpace) {
            continue;
        }

        if (spec->type == TRANS_MCFLIRT) {
            char *transPath = rsApplyTransformationGetMcFlirtTransformationPath(params->tmpDirPath, params->volumeIndex, spec->transformationId);
            rsStringAppend(callString, " -t ");
            rsStringAppend(callString, transPath);
            rsFree(transPath);
        } else if (spec->type == TRANS_ANTS) {
            rsStringAppend(callString, " -t ");
            rsStringAppend(callString, spec->file);
        } else if (spec->type == TRANS_FUGUE) {
            rsStringAppend(callString, " -t ");
            rsStringAppend(callString, spec->transOutPath);
        }
    }

    // execute
    if (params->verbose) {
        if (params->resourceTransformation) {
            fprintf(stdout, "Transforming resource file volume #%03d which will then be applied at the end of the transformation using:\n%s\n\n", (int) params->volumeIndex, callString);
        } else {
            fprintf(stdout, "Transforming final volume #%03d using:\n%s\n\n", (int) params->volumeIndex, callString);
        }
    }

    int returnStatus = system(callString);

    if (returnStatus != 0) {
        fprintf(stderr, "Error while converting volume %03d using:\n%s\n", (int)params->volumeIndex, callString);
        return FALSE;
    }

    rsFree(callString);

    return TRUE;
}

BOOL rsApplyTransformationApplyToVolume(const rsApplyTransformationApplyParams *params)
{
    // create tmp file for the input nifti
    char *inputName = rsMalloc(sizeof(char)*(strlen(params->tmpDirPath)+strlen(params->inputSuffix)+20));
    sprintf(inputName, "%s/%03d_%s.nii", params->tmpDirPath, params->volumeIndex, params->inputSuffix);
    char *outputName = rsMalloc(sizeof(char)*(strlen(params->tmpDirPath)+strlen(params->outputSuffix)+20));
    sprintf(outputName, "%s/%03d_%s.nii", params->tmpDirPath, params->volumeIndex, params->outputSuffix);
    char *inputMaskName = rsMalloc(sizeof(char)*(strlen(params->tmpDirPath)+strlen(params->inputMaskSuffix)+20));
    sprintf(inputMaskName, "%s/%03d_%s.nii", params->tmpDirPath, params->volumeIndex, params->inputMaskSuffix);
    char *outputMaskName = rsMalloc(sizeof(char)*(strlen(params->tmpDirPath)+strlen(params->outputMaskSuffix)+20));
    sprintf(outputMaskName, "%s/%03d_%s.nii", params->tmpDirPath, params->volumeIndex, params->outputMaskSuffix);
    char *paddedInputName = rsMalloc(sizeof(char)*(strlen(params->tmpDirPath)+strlen(params->paddedInputSuffix)+20));
    sprintf(paddedInputName, "%s/%03d_%s.nii", params->tmpDirPath, params->volumeIndex, params->paddedInputSuffix);
    rsNiftiFile *input = rsCloneNiftiFile(inputName, params->source, RSNIFTI_OPEN_ALLOC, 1);

    // write input nifti (single volume)
    double ***tmp = d3matrix(input->zDim-1, input->yDim-1, input->xDim-1);
    short indexToRead = params->volumeIndex;
    if (params->resourceTransformation && params->source->vDim == 1) {
        indexToRead = 0;
    }
    rsExtractVolumeFromRSNiftiFileBuffer(params->source, tmp[0][0], indexToRead);
    rsWriteVolumeToRSNiftiFileBuffer(input, tmp[0][0], 0);
    rsWriteNiftiHeader(input->fslio, "");
    FslWriteVolumes(input->fslio, input->data, input->vDim);
    rsCloseNiftiFileAndFree(input);
    rsFree(tmp[0][0]); rsFree(tmp[0]); rsFree(tmp);

    // pad input as parts of it will be removed by the lanczos-filter otherwise
    int nPaddingArguments = 10;
    char *paddingArguments[] = {
        "rszeropadding",
        rsStringConcat("--input=", inputName, NULL),
        rsStringConcat("--output=", paddedInputName, NULL),
        "--lx=5",
        "--ux=5",
        "--ly=5",
        "--uy=5",
        "--lz=5",
        "--uz=5",
        "--mirroredPadding"
    };
    rsZeropaddingParameters * paddingParams = rsZeropaddingParseParams(nPaddingArguments, paddingArguments);
    if (!paddingParams->parametersValid)
        return FALSE;
    rsZeropaddingInit(paddingParams);
    if (!paddingParams->parametersValid)
        return FALSE;
    rsZeropaddingRun(paddingParams);
    rsZeropaddingDestroy(paddingParams);

    // create a mask based on the padded input that will be used to remove the padding later
    rsNiftiFile *paddedInput = rsOpenNiftiFile(paddedInputName, RSNIFTI_OPEN_NONE);
    rsNiftiFile *inputMask = rsCloneNiftiFile(inputMaskName, paddedInput, RSNIFTI_OPEN_ALLOC, 1);
    tmp = d3matrix(inputMask->zDim-1, inputMask->yDim-1, inputMask->xDim-1);
    for (short x=0; x < inputMask->xDim; x++) {
        for (short y=0; y < inputMask->yDim; y++) {
            for (short z = 0; z < inputMask->zDim; z++) {
                if (x<5 || y<5 || z<5 || x>inputMask->xDim-6 || y>inputMask->yDim-6 || z>inputMask->zDim-6) {
                    tmp[z][y][x] = 0.0;
                } else {
                    tmp[z][y][x] = 1.0;
                }
            }
        }
    }
    rsWriteVolumeToRSNiftiFileBuffer(inputMask, tmp[0][0], 0);
    rsWriteNiftiHeader(inputMask->fslio, "");
    FslWriteVolumes(inputMask->fslio, inputMask->data, inputMask->vDim);
    rsCloseNiftiFileAndFree(inputMask);
    rsCloseNiftiFileAndFree(paddedInput);
    rsFree(tmp[0][0]); rsFree(tmp[0]); rsFree(tmp);

    // warp padded input
    if (!rsApplyTransformationRunANTs(params, paddedInputName, outputName, TRUE)) {
        return FALSE;
    }

    // warp padding mask
    if (!rsApplyTransformationRunANTs(params, inputMaskName, outputMaskName, FALSE)) {
        return FALSE;
    }

    // open resulting output nifti file and the warped mask
    rsNiftiFile *output = rsOpenNiftiFile(outputName, RSNIFTI_OPEN_READ);
    rsNiftiFile *outputMask = rsOpenNiftiFile(outputMaskName, RSNIFTI_OPEN_READ);

    if (!output->readable || !outputMask->readable) {
        fprintf(stderr, "Error: There was an error while executing antsApplyTransforms\n");
        return FALSE;
    }

    // read the single volume we got from ants
    double ***result = d3matrix(output->zDim-1, output->yDim-1, output->xDim-1);
    double ***resultMask = d3matrix(outputMask->zDim-1, outputMask->yDim-1, outputMask->xDim-1);
    rsExtractVolumeFromRSNiftiFileBuffer(output, result[0][0], 0);
    rsExtractVolumeFromRSNiftiFileBuffer(outputMask, resultMask[0][0], 0);

    // apply additional transformations that were already warped to the target space
    if (!params->resourceTransformation) {
        for (short i = 0; i < params->nTransformations; i++) {
            const rsApplyTransformationTransSpecification *spec = params->specs[i];

            if (!spec->isAppliedInTargetSpace) {
                continue;
            }

            if (spec->type == TRANS_MULTIPLICATION || spec->type == TRANS_DIVISION) {
                double ***tmp2 = d3matrix(output->zDim-1, output->yDim-1, output->xDim-1);
                rsExtractVolumeFromRSNiftiFileBuffer(spec->transOut, tmp2[0][0], params->volumeIndex);
                for (short x=0; x < output->xDim; x++) {
                    for (short y=0; y < output->yDim; y++) {
                        for (short z = 0; z < output->zDim; z++) {
                            if (spec->type == TRANS_MULTIPLICATION) {
                                result[z][y][x] *= tmp2[z][y][x];
                            } else {
                                result[z][y][x] /= tmp2[z][y][x];
                            }
                        }
                    }
                }
                rsFree(tmp2[0][0]); rsFree(tmp2[0]); rsFree(tmp2);
            }
        }
    }

    // apply padded mask to the output nifti
    for (short x=0; x < output->xDim; x++) {
        for (short y=0; y < output->yDim; y++) {
            for (short z = 0; z < output->zDim; z++) {
                if (resultMask[z][y][x] < 0.5) {
                    result[z][y][x] *= log(-1.0);
                }
            }
        }
    }

    // copy resulting volume to the appropriate index in the output volume
    rsWriteVolumeToRSNiftiFileBuffer(params->target, result[0][0], params->volumeIndex);
    rsCloseNiftiFileAndFree(output);
    rsCloseNiftiFileAndFree(outputMask);

    // cleanup
    if (!params->keepFiles) {
        unlink(inputName);
        unlink(outputName);
        unlink(paddedInputName);
        unlink(inputMaskName);
        unlink(outputMaskName);
    }

    rsFree(inputName);
    rsFree(outputName);
    rsFree(inputMaskName);
    rsFree(outputMaskName);
    rsFree(paddedInputName);
    rsFree(result[0][0]); rsFree(result[0]); rsFree(result);
    rsFree(resultMask[0][0]); rsFree(resultMask[0]); rsFree(resultMask);
}

void rsApplyTransformationDestroy(rsApplyTransformationParameters *p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }

    if ( p->output != NULL ) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }

    if ( p->transform != NULL ) {
        fclose(p->transform);
        p->transform = NULL;
    }

    rsApplyTransformationFreeParams(p);
}

char *rsApplyTransformationGetMcFlirtTransformationPath(const char *tmpDirPath, const int volumeIndex, const int transformationId)
{
    char *transPath = rsMalloc(sizeof(char)*(strlen(tmpDirPath) + 60));
    sprintf(transPath, "%s/%03d_%d_mcflirt_trans.txt", tmpDirPath, volumeIndex, transformationId);
    return transPath;
}

BOOL rsApplyTransformationConvertMcFlirtTransformations(rsNiftiFile* input, char ***transformations, size_t nVolumes, const char* transformationPath)
{
    // check if path is a directory
    struct stat s;
    BOOL isDirectory = stat(transformationPath, &s) == 0 && S_ISDIR(s.st_mode);

    if (!isDirectory) {
        fprintf(stderr, "The specified path to the MCFLIRT transformations does not point to a directory (path: %s)\n", transformationPath);
        return FALSE;
    }

    // load transformation files
    *transformations = (char**)rsMalloc(sizeof(char*)*nVolumes);
    for (short i=0; i<nVolumes; i++) {
        // assemble path
        char *path = (char*)rsMalloc(sizeof(char)*(strlen(transformationPath)+10));
        sprintf(path, "%s/MAT_%04d", transformationPath, (int)i);

        // load file
        double **mat;
        long nColumns, nRows;
        mat = rsLoadMatrixFromTxt(path, &nColumns, &nRows);

        if (mat == NULL) {
            return FALSE;
        }

        if (nColumns != 4 || nRows != 4) {
            fprintf(stderr, "Matrix format of file %s should have been 4x4, but was %dx%d!\n", path, nColumns, nRows);
            return FALSE;
        }

        // convert FSL mcflirt's transformation matrix to ANTs transformation format
        mat44 antsTrans; // <- ANTs transformation matrix
        mat44 M; // <- FSL's transformation matrix
        for (int i =0; i < 4; i++) {
            for (int j =0; j < 4; j++) {
                M.m[j][i] = mat[i][j];
            }
        }

        rsApplyTransformationConvertMcFlirtTransformMatrixToAntsTransformMatrix(&antsTrans, input, &M);

        // convert to ITK representation
        (*transformations)[i] = (char*)rsMalloc(sizeof(char*)*500);
        (*transformations)[i][0] = '\0';
        rsStringAppend((*transformations)[i], "#Insight Transform File V1.0\n");
        rsStringAppend((*transformations)[i], "#Transform 0\n");
        rsStringAppend((*transformations)[i], "Transform: AffineTransform_double_3_3\n");
        rsStringAppend((*transformations)[i], "Parameters:");
        for (short j=0; j<3; j++) {
            const size_t bufSize = strlen((*transformations)[i]);
            char *strBuf = &((*transformations)[i][bufSize]);
            sprintf(strBuf, " %.14f %.14f %.14f", antsTrans.m[j][0], antsTrans.m[j][1], antsTrans.m[j][2]);
        }
        for (short j=0; j<3; j++) {
            const size_t bufSize = strlen((*transformations)[i]);
            char *strBuf = &((*transformations)[i][bufSize]);
            sprintf(strBuf, " %.14f", antsTrans.m[j][3]);
        }
        rsStringAppend((*transformations)[i], "\n");
        {
            const size_t bufSize = strlen((*transformations)[i]);
            char *strBuf = &((*transformations)[i][bufSize]);
            sprintf(strBuf, "FixedParameters: %.14f %.14f %.14f\n", 0.0, 0.0, 0.0);
        }

        // cleanup
        rsFree(path);
        rsFree(mat[0]);
        rsFree(mat);
    }

    return TRUE;
}

void rsApplyTransformationConvertMcFlirtTransformMatrixToAntsTransformMatrix(mat44 *antsTransform, const rsNiftiFile* input, const mat44 *fslTransform)
{
    // declare some temporary variables
    gsl_matrix *tmp = gsl_matrix_alloc(4, 4);
    int signum;
    gsl_permutation *p;

    // M - FSL transformation matrix
    gsl_matrix *M = gsl_matrix_calloc(4, 4);
    for (int i =0; i < 4; i++) {
        for (int j =0; j < 4; j++) {
            gsl_matrix_set(M, i, j, fslTransform->m[i][j]);
        }
    }

    // R - the world coordinate matrix of the input nifti
    mat44 R44 = input->fslio->niftiptr->sform_code == NIFTI_XFORM_UNKNOWN
              ? input->fslio->niftiptr->qto_xyz
              : input->fslio->niftiptr->sto_xyz;
    gsl_matrix *R = gsl_matrix_calloc(4, 4);
    for (int i =0; i < 4; i++) {
        for (int j =0; j < 4; j++) {
            gsl_matrix_set(R, i, j, R44.m[i][j]);
        }
    }

    // S - a matrix containing the scaling factors of the input nifti
    gsl_matrix *S = gsl_matrix_calloc(4, 4);
    gsl_matrix_set(S, 0, 0, input->fslio->niftiptr->pixdim[1]);
    gsl_matrix_set(S, 1, 1, input->fslio->niftiptr->pixdim[2]);
    gsl_matrix_set(S, 2, 2, input->fslio->niftiptr->pixdim[3]);
    gsl_matrix_set(S, 3, 3, 1.0);

    // compute determinant of R
    p = gsl_permutation_calloc(4);
    gsl_matrix_memcpy(tmp , R);
    gsl_linalg_LU_decomp(tmp, p ,&signum);
    const double det = gsl_linalg_LU_det(tmp , signum);
    gsl_permutation_free(p);

    // if determinant > 0, flip x in S
    if (det > 0) {
        gsl_matrix_set(S, 0, 3, gsl_matrix_get(S, 0, 0) * (input->xDim-1));
        gsl_matrix_set(S, 0, 0, gsl_matrix_get(S, 0, 0) * -1.0);
    }

    // ras2lpi - RAS to LPI transformation (flip x and y)
    gsl_matrix *ras2lpi = gsl_matrix_alloc(4, 4); // -1  0  0  0
    gsl_matrix_set_identity(ras2lpi);             //  0 -1  0  0
    gsl_matrix_set(ras2lpi, 0, 0, -1.0);          //  0  0  1  0
    gsl_matrix_set(ras2lpi, 1, 1, -1.0);          //  0  0  0  1

    // compute inverse of S
    gsl_matrix *invS = gsl_matrix_alloc(4, 4);
    gsl_matrix_memcpy(tmp, S);
    p = gsl_permutation_alloc(4);
    gsl_linalg_LU_decomp(tmp, p, &signum);
    gsl_linalg_LU_invert(tmp, p, invS);
    gsl_permutation_free(p);

    // T = ras2lpi * R * inv(S)
    gsl_matrix *T = gsl_matrix_alloc(4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, R, invS, 0.0, tmp); // tmp <- R * inv(S)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ras2lpi, tmp, 0.0, T); // T = <- ras2lpi * tmp

    // compute inverse of T
    gsl_matrix *invT = gsl_matrix_alloc(4, 4);
    gsl_matrix_memcpy(tmp, T);
    p = gsl_permutation_alloc(4);
    gsl_linalg_LU_decomp(tmp, p, &signum);
    gsl_linalg_LU_invert(tmp, p, invT);
    gsl_permutation_free(p);

    // compute inverse of M
    gsl_matrix *invM = gsl_matrix_alloc(4, 4);
    gsl_matrix_memcpy(tmp, M);
    p = gsl_permutation_alloc(4);
    gsl_linalg_LU_decomp(tmp, p, &signum);
    gsl_linalg_LU_invert(tmp, p, invM);
    gsl_permutation_free(p);

    // antsTrans = T * M * inv(T)
    gsl_matrix *antsTrans = gsl_matrix_alloc(4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invM, invT, 0.0, tmp); // tmp <- M * inv(T)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, tmp, 0.0, antsTrans); // antsTrans <- T * tmp

    // copy result
    for (int i =0; i < 4; i++) {
        for (int j =0; j < 4; j++) {
            antsTransform->m[i][j] = gsl_matrix_get(antsTrans, i, j);
        }
    }

    // cleanup
    gsl_matrix_free(tmp);
    gsl_matrix_free(M);
    gsl_matrix_free(invM);
    gsl_matrix_free(R);
    gsl_matrix_free(S);
    gsl_matrix_free(invS);
    gsl_matrix_free(T);
    gsl_matrix_free(invT);
    gsl_matrix_free(ras2lpi);
    gsl_matrix_free(antsTrans);
}

BOOL rsApplyTransformationConvertFugueShiftToANTsWarp(const rsNiftiFile* input, const rsNiftiFile* shift, const char* warpPath, BOOL verbose)
{
    rsNiftiFile *warp = rsCloneNiftiFile(warpPath, shift, RSNIFTI_OPEN_ALLOC, 3);
    nifti_image *warpImage = warp->fslio->niftiptr;

    // Determine phase encoding direction
    rsNiftiExtendedHeaderInformation* info;
    info = rsNiftiFindExtendedHeaderInformation(input->fslio->niftiptr);

    if (info == NULL || info->PhaseEncodingDirection == NULL || strlen(&info->PhaseEncodingDirection[0]) != 2) {
        fprintf(stderr, "Could not determine phase encoding direction from the nifti header of %s\n", input->path);
        return FALSE;
    }

    char *phaseEncDir = (char*)rsMalloc(sizeof(char)*3);
    sprintf(phaseEncDir, "%s", info->PhaseEncodingDirection);

    if (phaseEncDir[0]!='x' && phaseEncDir[0]!='y' && phaseEncDir[0]!='z') {
        fprintf(stderr, "The dimension of the phase encoding direction could not be read from the nifti header of %s\n", input->path);
        return FALSE;
    }

    if (phaseEncDir[1]!='+' && phaseEncDir[1]!='-') {
        fprintf(stderr, "The sign of the phase encoding direction could not be read from the nifti header of %s\n", input->path);
        return FALSE;
    }

    short phaseEncDim = 0;
    if (phaseEncDir[0]=='y') {
        phaseEncDim = 1;
    } else if (phaseEncDir[0]=='z') {
        phaseEncDim = 2;
    }

    double phaseEncSign = phaseEncDir[1] == '-' ? -1.0 : +1.0;

    // Convert world matrix from RAS to LPI
    mat44 worldMatrix = warpImage->sform_code == NIFTI_XFORM_UNKNOWN
                        ? warpImage->qto_xyz
                        : warpImage->sto_xyz;

    // R - the world coordinate matrix of the input nifti
    gsl_matrix *R = gsl_matrix_calloc(4, 4);
    for (int i =0; i < 4; i++) {
        for (int j =0; j < 4; j++) {
            gsl_matrix_set(R, i, j, worldMatrix.m[i][j]);
        }
    }

    // ras2lpi - RAS to LPI transformation (flip x and y)
    gsl_matrix *ras2lpi = gsl_matrix_alloc(4, 4); // -1  0  0  0
    gsl_matrix_set_identity(ras2lpi);             //  0 -1  0  0
    gsl_matrix_set(ras2lpi, 0, 0, -1.0);          //  0  0  1  0
    gsl_matrix_set(ras2lpi, 1, 1, -1.0);          //  0  0  0  1

    // antsR = ras2lpi * R
    gsl_matrix *antsR = gsl_matrix_alloc(4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ras2lpi, R, 0.0, antsR); // antsR = <- ras2lpi * R
    gsl_matrix_free(R);
    gsl_matrix_free(ras2lpi);

    if (verbose) {
        fprintf(stdout, "Phase encoding direction: %s\n", phaseEncDir);
        fprintf(stdout, "Multiplying shift with: %.0f\n", phaseEncSign);
        fprintf(stdout, "World matrix of the supplied fugue voxel shift:\n");
        for (short i = 0; i < 4; i++) {
            fprintf(stdout, " %+.6f  %+.6f  %+.6f  %+.6f\n", worldMatrix.m[i][0], worldMatrix.m[i][1], worldMatrix.m[i][2], worldMatrix.m[i][3]);
        }
        fprintf(stdout, "World matrix used for converting %s shifts to ANTs warp vectors:\n", phaseEncDir);
        for (short i = 0; i < 4; i++) {
            fprintf(stdout, " %+.6f  %+.6f  %+.6f  %+.6f\n", gsl_matrix_get(antsR, i, 0), gsl_matrix_get(antsR, i, 1), gsl_matrix_get(antsR, i, 2), gsl_matrix_get(antsR, i, 3));
        }
    }

    // Create ANTs warp from FSL's voxel shift map
    double ***shiftData = d3matrix(shift->zDim-1, shift->yDim-1, shift->xDim-1);
    double ***warpData = d3matrix(shift->zDim-1, shift->yDim-1, shift->xDim-1);
    rsExtractVolumeFromRSNiftiFileBuffer(shift, shiftData[0][0], 0);
    for (short dim=0; dim<3; dim++) {
        for (short z = 0; z < shift->zDim; z++) {
            for (short y=0; y < shift->yDim; y++) {
                for (short x=0; x < shift->xDim; x++) {

                    const double phaseCorrectedShift = shiftData[z][y][x] * phaseEncSign;

                    const double shiftVX[3] = {
                        (phaseEncDim == 0) ? phaseCorrectedShift : 0.0,
                        (phaseEncDim == 1) ? phaseCorrectedShift : 0.0,
                        (phaseEncDim == 2) ? phaseCorrectedShift : 0.0
                    };

                    double shiftMM[3] = {0.0, 0.0, 0.0};
                    for (short i=0; i<3; i++) {
                        for (short j=0; j<3; j++) {
                            shiftMM[i] += gsl_matrix_get(antsR, i, j) * shiftVX[j];
                        }
                    }

                    warpData[z][y][x] = shiftMM[dim];
                }
            }
        }
        rsWriteVolumeToRSNiftiFileBuffer(warp, warpData[0][0], dim);
    }
    rsFree(warpData[0][0]); rsFree(warpData[0]); rsFree(warpData);
    rsFree(shiftData[0][0]); rsFree(shiftData[0]); rsFree(shiftData);
    gsl_matrix_free(antsR);

    // Modify dimensionality of the header so that the 4th dimension moves to the 5th
    FslSetDimensionality(warp->fslio, 5);
    warpImage->nt = 1; // 4th dim
    warpImage->nu = 3; // 5th dim
    warpImage->nv = 1; // 6th dim
    warpImage->nw = 1; // 7th dim

    warpImage->dim[0] = 5;
    warpImage->dim[1] = warpImage->nx;
    warpImage->dim[2] = warpImage->ny;
    warpImage->dim[3] = warpImage->nz;
    warpImage->dim[4] = warpImage->nt;
    warpImage->dim[5] = warpImage->nu;
    warpImage->dim[6] = warpImage->nv;
    warpImage->dim[7] = warpImage->nw;

    FslSetIntent(warp->fslio, NIFTI_INTENT_VECTOR, 0.0f, 0.0f, 0.0f);

    // Write out warp file
    rsWriteNiftiHeader(warp->fslio, "");
    FslWriteVolumes(warp->fslio, warp->data, warp->vDim);
    rsCloseNiftiFileAndFree(warp);

    return TRUE;
}

void rsApplyTransformationLoadTransformationFile(char ***output, size_t *nTransformations, FILE *file)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    *nTransformations = 0;

    // determine number of transformations
    while ((read = getline(&line, &len, file)) != -1) {
        *nTransformations = *nTransformations + 1;
    }

    *output = (char**)rsMalloc(sizeof(char*) * *nTransformations);
    fseek(file, 0, SEEK_SET);

    // determine length of indidvidual transformations
    short i = 0;
    while ((read = getline(&line, &len, file)) != -1) {
        (*output)[i] = (char*)rsMalloc(sizeof(char)*(len+1));
        i++;
    }

    fseek(file, 0, SEEK_SET);

    // read in transformations
    i = 0;
    while ((read = getline(&line, &len, file)) != -1) {
        sprintf((*output)[i], "%s", line);
        const size_t lastCharPos = strlen(line)-1;
        const char lastChar = (*output)[i][lastCharPos];
        if (lastChar == '\n' || lastChar == '\r') {
            (*output)[i][lastCharPos] = '\0';
        }
        i++;
    }
}

BOOL rsApplyTransformationParseTransformationFile(rsApplyTransformationTransSpecification*** transformations, char **transformationFile, const size_t nTransformations)
{
    *transformations = (rsApplyTransformationTransSpecification**)rsMalloc(sizeof(rsApplyTransformationTransSpecification*)*nTransformations);

    for (size_t i=0; i<nTransformations; i++) {
        (*transformations)[i] = (rsApplyTransformationTransSpecification*)rsMalloc(sizeof(rsApplyTransformationTransSpecification));
        size_t preamble;

        if (rsStringStartsWith(transformationFile[i], "-mcflirt ")) {
            (*transformations)[i]->type = TRANS_MCFLIRT;
            preamble = strlen("-mcflirt ");
            (*transformations)[i]->isAppliedInTargetSpace = FALSE;
        } else if (rsStringStartsWith(transformationFile[i], "-ants ")) {
            (*transformations)[i]->type = TRANS_ANTS;
            (*transformations)[i]->isAppliedInTargetSpace = FALSE;
            preamble = strlen("-ants ");
        } else if (rsStringStartsWith(transformationFile[i], "-mult ")) {
            (*transformations)[i]->type = TRANS_MULTIPLICATION;
            (*transformations)[i]->isAppliedInTargetSpace = TRUE;
            preamble = strlen("-mult ");
        } else if (rsStringStartsWith(transformationFile[i], "-div ")) {
            (*transformations)[i]->type = TRANS_DIVISION;
            (*transformations)[i]->isAppliedInTargetSpace = TRUE;
            preamble = strlen("-div ");
        } else if (rsStringStartsWith(transformationFile[i], "-fugue ")) {
            (*transformations)[i]->type = TRANS_FUGUE;
            (*transformations)[i]->isAppliedInTargetSpace = FALSE;
            preamble = strlen("-fugue ");
        } else {
            fprintf(stderr, "Could not parse the following transformation: \"%s\"\n", transformationFile[i]);
            return FALSE;
        }
        (*transformations)[i]->transformationId = i;
        const size_t fileNameSize = strlen(&transformationFile[i][preamble]);
        (*transformations)[i]->file = (char*)rsMalloc(sizeof(char)*(fileNameSize+1));
        sprintf((*transformations)[i]->file, "%s", &transformationFile[i][preamble]);
    }

    return TRUE;
}

char *rsApplyTransformationGetANTsPath()
{
    const char *confPath = CONFIG_PATH"/rstools.conf";
    FILE *config = fopen(confPath, "r");

    if (config == NULL) {
        fprintf(stderr, "Could not read rstools config file to retrieve the path to ANTs (%s)\n", confPath);
        return NULL;
    }

    int length;
    char *line = rsMalloc(sizeof(char)*1000);
    const char* pathPreamble = "ANTSPATH=";
    char *antsPath = NULL;

    while (rsReadline(config, line, &length)) {
        if (!rsStringStartsWith(line, pathPreamble)) {
            continue;
        }
        antsPath = rsString(&line[strlen(pathPreamble)]);
    }

    fclose(config);
    rsFree(line);

    if (antsPath == NULL) {
        fprintf(stderr, "The rstools config file did not contain the path to ANTs (%s)\n", confPath);
        return NULL;
    }

    return antsPath;
}