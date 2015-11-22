#include "rsdeoblique_ui.h"

rsDeobliqueParameters *rsDeobliqueInitParameters() {
    rsDeobliqueParameters *p = (rsDeobliqueParameters*)rsMalloc(sizeof(rsDeobliqueParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->transformationpath   = NULL;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->transform            = NULL;
    p->parametersValid      = FALSE;
    
    return p;
}

void rsDeobliqueFreeParams(rsDeobliqueParameters *p) {
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->transformationpath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->transform);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsDeobliqueParameters *rsDeobliqueParseParams(int argc, char * argv[]) {

    rsDeobliqueParameters *p = rsDeobliqueInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsDeobliqueBuildInterface(p);
    
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

    p->parametersValid = parsingSuccessful;
    return p;
}

void rsDeobliqueBuildInterface(rsDeobliqueParameters *p)
{  
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Resamples the volume such that ";

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
    o->name                = "retainVoxelSizes";
    o->storage             = &p->retainVoxelSizes;
    o->cli_description     = "only orthogonalize matrix, but switch to an isotropic resolution";
    rsUIAddOption(p->interface, o);
	
    o = rsUINewOption();
    o->name                = "transformation";
    o->shorthand           = 'r';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->transformationpath;
    o->cli_description     = "if specified, an ITK transformation file will be saved to this file path describing the input to output transformation. This file can be used with ANTs.";
    o->cli_arg_description = "<*.txt>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
}
