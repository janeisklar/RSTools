#include "rssmoothing_ui.h"

rsSmoothingParameters *rsSmoothingInitParameters() {
    rsSmoothingParameters *p = (rsSmoothingParameters*)rsMalloc(sizeof(rsSmoothingParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->kernelSizeFWHM       = -1.0;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->threads              = 1;
    p->progressCallback     = NULL;
    p->parametersValid      = FALSE;
    
    return p;
}

void rsSmoothingFreeParams(rsSmoothingParameters *p) {
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsSmoothingParameters *rsSmoothingParseParams(int argc, char * argv[]) {

    rsSmoothingParameters *p = rsSmoothingInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsSmoothingBuildInterface(p);
    
    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    // check if the required arguments have been provided
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified!\n");
        return p;
    }
    
    if ( p->outputpath == NULL ) {
        fprintf(stderr, "An output volume must be specified!\n");
        return p;
    }
    
    if ( p->kernelSizeFWHM < 0.000001 ) {
        fprintf(stderr, "A kernel size must be specified!\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsSmoothingBuildInterface(rsSmoothingParameters *p)
{  
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Applies gaussian smoothing for a given kernel size to a 4D nifti.";

    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the 4D volume to be smoothed";
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
    o->name                = "kernelSize";
    o->shorthand           = 'k';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->kernelSizeFWHM;
    o->cli_description     = "The size (sigma) of the gaussian kernel in mm (FWHM)";
    o->cli_arg_description = "<float>";
    o->defaultValue        = "6.0";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "threads";
    o->shorthand           = 't';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->threads;
    o->cli_description     = "number of threads used for processing";
    o->cli_arg_description = "<n>";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
	
    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
}
