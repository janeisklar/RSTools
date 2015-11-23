#ifndef rstools_applytransformation_ui_h
#define rstools_applytransformation_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

enum rsApplyTransformationTransType {
    TRANS_MCFLIRT,
    TRANS_ANTS,
    TRANS_MULTIPLICATION,
    TRANS_DIVISION,
    TRANS_FUGUE
};

typedef struct {
    enum rsApplyTransformationTransType type;
    char* file;
    char** mcflirtAntsFiles;
    rsNiftiFile *transIn;
    rsNiftiFile *transOut;
    char *transOutPath;
    BOOL isAppliedInTargetSpace;
    short transformationId;
} rsApplyTransformationTransSpecification;

typedef struct {
    char *inputpath;
    char *outputpath;
    char *transformationpath;
    char *referencepath;
    char *headerReferencePath;
    char *antsPath;

    char *callString;

    int threads;
    BOOL verbose;
    BOOL keepFiles;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *output;
    FILE *transform;

    rsApplyTransformationTransSpecification** specs;
    size_t nTransformations;

    rsUIInterface *interface;

} rsApplyTransformationParameters;

rsApplyTransformationParameters *rsApplyTransformationParseParams(int argc, char * argv[]);
rsApplyTransformationParameters *rsApplyTransformationInitParameters();
void rsApplyTransformationBuildInterface(rsApplyTransformationParameters *p);
void rsApplyTransformationFreeParams(rsApplyTransformationParameters *p);

#ifdef __cplusplus
}
#endif

#endif
