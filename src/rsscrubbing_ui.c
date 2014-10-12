#include "rsscrubbing_ui.h"

rsScrubbingParameters *rsScrubbingInitParameters() {
    rsScrubbingParameters *p = (rsScrubbingParameters*)rsMalloc(sizeof(rsScrubbingParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->flaggedpath          = NULL;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->parametersValid      = FALSE;
    
    return p;
}

void rsScrubbingFreeParams(rsScrubbingParameters *p) {
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->flaggedpath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsScrubbingParameters *rsScrubbingParseParams(int argc, char * argv[]) {

    rsScrubbingParameters *p = rsScrubbingInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsScrubbingBuildInterface(p);
    
    // parse
    p->parametersValid = rsUIParse(p->interface, argc, argv, (void*)p);
    
    return p;
}

void rsScrubbingBuildInterface(rsScrubbingParameters *p)
{  
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Given a 4D-nifti and a txt-file(such as --flagged from rsmotionscrubbing) containing indices referring to volumes that are to be removed from the 4D-nifti.";

    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the 4D volume to be scrubbed";
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
    o->name                = "flagged";
    o->shorthand           = 'e';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->flaggedpath;
    o->cli_description     = "file containing the indices of all flagged volumes that will be removed from the input (indices should be 0-based and in the rows)";
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
