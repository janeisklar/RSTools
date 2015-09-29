#include "rsinfo_ui.h"

rsInfoParameters* rsInfoInitParameters()
{
    rsInfoParameters *p = (rsInfoParameters*)rsMalloc(sizeof(rsInfoParameters));
    
    p->inputpath            = NULL;
    p->dicompath            = NULL;
    p->parametersValid      = FALSE;
    p->input                = NULL;
    p->dicom                = NULL;
    p->showComments         = FALSE;
    p->showInfo             = FALSE;
    p->infoKey              = NULL;
    p->interface            = NULL;

    return p;
}

rsInfoParameters* rsInfoParseParams(int argc, char * argv[])
{
    rsInfoParameters *p = rsInfoInitParameters();
    rsInfoBuildInterface(p);
    
    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    // check if the required arguments have been provided
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified(--input)!\n");
        return p;
    }

    p->parametersValid = parsingSuccessful;
    return p;
}

void rsInfoBuildInterface(rsInfoParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Given a Nifti that has been created by any of the rstools, this tool reads out the header information it entails.";
    
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the input volume to be regressed";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "dicom";
    o->shorthand           = 'd';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->dicompath;
    o->cli_description     = "save the dicom in the nifti header to the specified location";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "key";
    o->shorthand           = 'k';
    o->type                = G_OPTION_ARG_STRING;
    o->storage             = &p->infoKey;
    o->cli_description     = "key that will be read out from the extended nifti header information (overwrites --info). For a list of keys see --info";
    o->cli_arg_description = "<key>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "info";
    o->shorthand           = 'n';
    o->storage             = &p->showInfo;
    o->cli_description     = "show extended nifti header information";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "comments";
    o->shorthand           = 'c';
    o->storage             = &p->showComments;
    o->cli_description     = "show the comments in the nifti header";
    rsUIAddOption(p->interface, o);
}

void rsInfoFreeParams(rsInfoParameters* p)
{
    rsFree(p->inputpath);
    rsFree(p->dicompath);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}
