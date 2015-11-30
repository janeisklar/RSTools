#include <src/utils/rsui.h>
#include "rsapplytransformation_ui.h"

rsApplyTransformationParameters *rsApplyTransformationInitParameters() {
    rsApplyTransformationParameters *p = (rsApplyTransformationParameters*)rsMalloc(sizeof(rsApplyTransformationParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->referencepath        = NULL;
    p->headerReferencePath  = NULL;
    p->antsPath             = NULL;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->keepFiles            = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->threads              = 1;
    p->transform            = NULL;
    p->nTransformations     = 0;
    p->specs                = NULL;
    p->parametersValid      = FALSE;
    p->progressCallback     = NULL;
    
    return p;
}

void rsApplyTransformationFreeParams(rsApplyTransformationParameters *p) {
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->referencepath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->transform);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsApplyTransformationParameters *rsApplyTransformationParseParams(int argc, char * argv[]) {

    rsApplyTransformationParameters *p = rsApplyTransformationInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsApplyTransformationBuildInterface(p);
    
    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    // check if the required arguments have been provided
    if (p->inputpath == NULL) {
        fprintf(stderr, "No input volume specified!\n");
        return p;
    }
    
    if (p->outputpath == NULL) {
        fprintf(stderr, "An output volume must be specified!\n");
        return p;
    }

    if (p->referencepath == NULL) {
        fprintf(stderr, "A reference volume must be specified!\n");
        return p;
    }

    if (p->transformations == NULL || g_strv_length(p->transformations) < 1) {
        fprintf(stderr, "At least one transformation needs to be specified!\n");
        return p;
    }
    p->nTransformations = g_strv_length(p->transformations);

    p->parametersValid = parsingSuccessful;
    return p;
}

void rsApplyTransformationBuildInterface(rsApplyTransformationParameters *p)
{  
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Applies a complex sequence of various transformations using ANTs as speecified by the --trans/-T parameter.";

    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the 4D volume to processed";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "output";
    o->shorthand           = 'o';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->outputpath;
    o->cli_description     = "the volume in which the result will be saved";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "reference";
    o->shorthand           = 'r';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->referencepath;
    o->cli_description     = "a reference volume that has the desired orientation/scaling";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "headerReference";
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->headerReferencePath;
    o->cli_description     = "a volume from which the extended header information will be copied (world matrices, voxel sizes and dim lengths will be taken from the volume supplied with --reference/-r)";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "trans";
    o->shorthand           = 'T';
    o->type                = G_OPTION_ARG_STRING_ARRAY;
    o->nLines              = 10;
    o->storage             = &p->transformations;
    o->cli_description     = "Specifies one transformation which is applied using antsApplyTransforms in the reversed order of occurrence, therefore the first transformation that is specified is also the first one that is applied to the data (which corresponds to the last transformation in ANT's notation). Transformations should be given in the following format \"-T type,path\" or \"--trans=type,path\", where 'type' is the type of the transformation (i.e. ants, mcflirt, fugue, mult, div) and 'path' the path to a transformation file or folder. In the case of 'mcflirt', the 'path' is assumed to be a folder with rotation matrices. The types 'mult' and 'div' can be used to specify a multiplication/division with a constant 3D or 4D nifti such as the estimated bias field, which will be warped as well and applied after interpolating the input nifti in the same space. The transformation type 'fugue' can be used to specify a voxel shift map for distortion correction.";
    o->gui_description     = "Specify one transformation per line which will be applied using antsApplyTransforms in the reversed order of occurrence, therefore the first transformation that is specified is also the first one that is applied to the data (which corresponds to the last transformation in ANT's notation). Transformations should be given in the following format \"type,path\", where 'type' is the type of the transformation (i.e. ants, mcflirt, fugue, mult, div) and 'path' the path to a transformation file or folder. In the case of 'mcflirt', the 'path' is assumed to be a folder with rotation matrices. The types 'mult' and 'div' can be used to specify a multiplication/division with a constant 3D or 4D nifti such as the estimated bias field, which will be warped as well and applied after interpolating the input nifti in the same space. The transformation type 'fugue' can be used to specify a voxel shift map for distortion correction.";
    o->cli_arg_description = "<type,param>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "threads";
    o->shorthand           = 't';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->threads;
    o->cli_description     = "number of threads used for processing. a value greater one disables some of the debug information that were activated using (-verbose)";
    o->cli_arg_description = "<n>";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "keepFiles";
    o->shorthand           = 'k';
    o->storage             = &p->keepFiles;
    o->cli_description     = "(DEBUG) does not delete temporary files that were created while applying the transformations. (use together with --verbose)";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
}
