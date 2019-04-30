#include "rsinfo_ui.h"

rsInfoParameters* rsInfoInitParameters()
{
    rsInfoParameters *p = (rsInfoParameters*)rsMalloc(sizeof(rsInfoParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->dicompath            = NULL;
    p->extensionSource      = NULL;
    p->parametersValid      = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->extensionInput       = NULL;
    p->dicom                = NULL;
    p->showComments         = FALSE;
    p->showInfo             = FALSE;
    p->infoKey              = NULL;
    p->modArgs              = NULL;
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
    o->cli_description     = "the input volume from which the header information will be read out";
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
    o->cli_description     = "key that will be read out  from the extended nifti header information (overwrites --info). For a list of keys see --info. The keys are not case-sensitive.";
    o->cli_arg_description = "<key>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "m";
    o->shorthand           = 'm';
    o->type                = G_OPTION_ARG_STRING_ARRAY;
    o->storage             = &p->modArgs;
    o->cli_description     = "key/value pairs that define header properties that will be modified in-place (--input / -i). The keys are not case-sensitive. Example: '-i foo.nii -m PatientAge=045Y' changes the patient age in the header of foo.nii. It can be specified more than once to modify multiple header properties. To set a string containing spaces, use quotes around the value (f.e. -m InstitutionName=\"My Institution Name\"). Arrays should be supplied as a list of comma-separated values (f.e. -m MosaicRefAcqTimes=1.23,4.56,7.89,...)";
    o->cli_arg_description = "<key=val>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "copyExtension";
    o->shorthand           = 'x';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->extensionSource;
    o->cli_description     = "allows to specify a nifti from which all extensions will be copied over to the output. (requires the specification of --output / -o)";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "output";
    o->shorthand           = 'o';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->outputpath;
    o->cli_description     = "instead of applying the changes in-place, a modified copy of the whole file will be created at the specified path";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "comments";
    o->shorthand           = 'c';
    o->storage             = &p->showComments;
    o->cli_description     = "show only the comments in the nifti header";
    rsUIAddOption(p->interface, o);
}

void rsInfoFreeParams(rsInfoParameters* p)
{
    rsFree(p->inputpath);
    rsFree(p->dicompath);
    rsFree(p->extensionSource);
    rsFree(p->outputpath);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}
