#include "rsrescale_ui.h"

rsRescaleParameters *rsRescaleInitParameters() {
    rsRescaleParameters *p = (rsRescaleParameters*)rsMalloc(sizeof(rsRescaleParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->scale[0]             = -1.0;
    p->scale[1]             = -1.0;
    p->scale[2]             = -1.0;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->threads              = 1;
    p->progressCallback     = NULL;
    p->parametersValid      = FALSE;
    
    return p;
}

void rsRescaleFreeParams(rsRescaleParameters *p) {
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsRescaleParameters *rsRescaleParseParams(int argc, char * argv[]) {

    rsRescaleParameters *p = rsRescaleInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsRescaleBuildInterface(p);
    
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
    
    if ( p->scale[0] < 0.000001 ) {
        fprintf(stderr, "A scaling factor for the voxel dimensions in x-direction needs to be specified!!\n");
        return p;
    }

    if ( p->scale[1] < 0.000001 ) {
        fprintf(stderr, "A scaling factor for the voxel dimensions in y-direction needs to be specified!!\n");
        return p;
    }

    if ( p->scale[2] < 0.000001 ) {
        fprintf(stderr, "A scaling factor for the voxel dimensions in z-direction needs to be specified!!\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsRescaleBuildInterface(rsRescaleParameters *p)
{  
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Spatially rescales the given nifti using a Lanczos interpolation filter";

    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the 4D volume to be rescaled";
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
    o->name                = "scaleX";
    o->shorthand           = 'x';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->scale[0];
    o->cli_description     = "the factor by which to scale the voxel dimensions in x-directiom";
    o->cli_arg_description = "<float>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "scaleY";
    o->shorthand           = 'y';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->scale[1];
    o->cli_description     = "the factor by which to scale the voxel dimensions in y-directiom";
    o->cli_arg_description = "<float>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "scaleZ";
    o->shorthand           = 'z';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->scale[2];
    o->cli_description     = "the factor by which to scale the voxel dimensions in z-directiom";
    o->cli_arg_description = "<float>";
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
