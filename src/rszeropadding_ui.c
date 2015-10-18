#include "rszeropadding_ui.h"

rsZeropaddingParameters *rsZeropaddingInitParameters() {
    rsZeropaddingParameters *p = (rsZeropaddingParameters*)rsMalloc(sizeof(rsZeropaddingParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->padding[0]           = 0;
    p->padding[1]           = 0;
    p->padding[2]           = 0;
    p->padding[3]           = 0;
    p->padding[4]           = 0;
    p->padding[5]           = 0;
    p->paddingValue         = 0.0;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->parametersValid      = FALSE;
    
    return p;
}

void rsZeropaddingFreeParams(rsZeropaddingParameters *p) {
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsZeropaddingParameters *rsZeropaddingParseParams(int argc, char * argv[]) {

    rsZeropaddingParameters *p = rsZeropaddingInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsZeropaddingBuildInterface(p);
    
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

void rsZeropaddingBuildInterface(rsZeropaddingParameters *p)
{  
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Allows to pad the volume by a specified number of voxels on the corresponding side. The value that the newly created space is filled with can be specified. It may also be used for cropping by supplying a negative padding. The world matrices (sform/qform) are adjusted accordingly.";

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
    o->name                = "lx";
    o->shorthand           = 'a';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->padding[0];
    o->cli_description     = "number of voxels the lower end of x is padded";
    o->cli_arg_description = "<int>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "ux";
    o->shorthand           = 'b';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->padding[1];
    o->cli_description     = "number of voxels the upper end of x is padded";
    o->cli_arg_description = "<int>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "ly";
    o->shorthand           = 'c';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->padding[2];
    o->cli_description     = "number of voxels the lower end of y is padded";
    o->cli_arg_description = "<int>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "uy";
    o->shorthand           = 'd';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->padding[3];
    o->cli_description     = "number of voxels the upper end of y is padded";
    o->cli_arg_description = "<int>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "lz";
    o->shorthand           = 'e';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->padding[4];
    o->cli_description     = "number of voxels the lower end of z is padded";
    o->cli_arg_description = "<int>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "uz";
    o->shorthand           = 'f';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->padding[5];
    o->cli_description     = "number of voxels the upper end of z is padded";
    o->cli_arg_description = "<int>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "value";
    o->shorthand           = 'p';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->paddingValue;
    o->cli_description     = "the value that padded is used for the padded space";
    o->cli_arg_description = "<double>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
}
