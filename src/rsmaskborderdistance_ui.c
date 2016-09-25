#include "rsmaskborderdistance_ui.h"

rsMaskBorderDistanceParameters *rsMaskBorderDistanceInitParameters()
{
    rsMaskBorderDistanceParameters *p = (rsMaskBorderDistanceParameters*)rsMalloc(sizeof(rsMaskBorderDistanceParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->parametersValid      = FALSE;
    p->threads              = 1;
    p->interface            = NULL;
    p->progressCallback     = NULL;
    
    return p;
}

void rsMaskBorderDistanceFreeParams(rsMaskBorderDistanceParameters *p)
{
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsMaskBorderDistanceParameters *rsMaskBorderDistanceParseParams(int argc, char * argv[])
{

    rsMaskBorderDistanceParameters *p = rsMaskBorderDistanceInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsMaskBorderDistanceBuildInterface(p);
    
    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified(--input)!\n");
        return p;
    }
    
    if ( p->outputpath == NULL ) {
        fprintf(stderr, "An output path for the filtered data must be specified(--filtered)!\n");
        return p;
    }

    p->parametersValid = parsingSuccessful;
    return p;
}

void rsMaskBorderDistanceBuildInterface(rsMaskBorderDistanceParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Takes in a mask and replaces every voxel in the mask (>0) with the computed distance (in mm) to the nearest non-mask voxel.";
    
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the input volume that is going to be filtered";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "output";
    o->shorthand           = 'o';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->outputpath;
    o->cli_description     = "the output volume in which the filtered data will be saved";
    o->cli_arg_description = "<volume>";
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
