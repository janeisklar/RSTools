#include "rsmotionscrubbing_ui.h"

rsMotionScrubbingParameters *rsMotionScrubbingInitParameters() {
    rsMotionScrubbingParameters *p = (rsMotionScrubbingParameters*)rsMalloc(sizeof(rsMotionScrubbingParameters));
    
    p->inputpath            = NULL;
    p->maskpath             = NULL;
    p->outputpath           = NULL;
    p->realignmentpath      = NULL;
    p->fdpath               = NULL;
    p->dvarspath            = NULL;
    p->flaggedpath          = NULL;
    p->callString           = NULL;
    p->fdthreshold          = 1.0;
    p->dvarsthreshold       = 0.05;
    p->useModal             = FALSE;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->parametersValid      = FALSE;
    p->mask                 = NULL;
    p->interface            = NULL;
    p->rp                   = NULL;
    p->rpformat             = RSTOOLS_REALIGNMENT_PARAMETER_FORMAT_SPM;
    p->maskPoints           = NULL;
    
    return p;
}

void rsMotionScrubbingFreeParams(rsMotionScrubbingParameters *p) {
    rsFree(p->inputpath);
    rsFree(p->maskpath);
    rsFree(p->outputpath);
    rsFree(p->realignmentpath);
    rsFree(p->fdpath);
    rsFree(p->dvarspath);
    rsFree(p->flaggedpath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->callString);
    rsFree(p->mask);
    if ( p->rp != NULL ) {
        rsFree(p->rp[0]);
    }
    rsFree(p->rp);
    rsFree(p->maskPoints);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsMotionScrubbingParameters *rsMotionScrubbingParseParams(int argc, char * argv[]) {

    rsMotionScrubbingParameters *p = rsMotionScrubbingInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsMotionScrubbingBuildInterface(p);
    
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    // check if the required arguments have been provided
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified!\n");
        return p;
    }
    
    if ( p->maskpath == NULL ) {
        fprintf(stderr, "A binary mask must be specified!\n");
        return p;
    }
    
    if ( p->realignmentpath == NULL ) {
        fprintf(stderr, "A realignment parameter file must be specified!\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsMotionScrubbingBuildInterface(rsMotionScrubbingParameters *p)
{  
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Given a 4D-Nifti and a txt-file containing the 6 head realignment parameters, this application performs motion-scrubbing as described in: Power, Jonathan D., et al. \"Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion.\" Neuroimage 59.3 (2012): 2142-2154. APA";

    GOptionArgFunc cbRPFormat = (GOptionArgFunc)rsMotionScrubbingParseRPFormat;

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
    o->name                = "rp";
    o->shorthand           = 'r';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->realignmentpath;
    o->cli_description     = "the file containing the realignment parameters";
    o->cli_arg_description = "<*.txt>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "rpformat";
    o->shorthand           = 'f';
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbRPFormat;
    o->cli_description     = "the format of the supplied realignment parameters (defaults to SPM)";
    o->cli_arg_description = "<format>";
    o->defaultValue        = "spm";
    rsUIOptionValue allowedValues[] = {
      {"spm",  "Column order: x,y,z,pitch,roll,yaw. Rotations are supplied in degrees."},
      {"fsl",  "Column order: pitch,roll,yaw,x,y,z. Rotations are supplied in radians."},
      NULL
    };
    rsUISetOptionValues(o, allowedValues);
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "dvarsthreshold";
    o->shorthand           = 'a';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->dvarsthreshold;
    o->cli_description     = "(optional) DVARs threshold";
    o->cli_arg_description = "<float>";
    o->defaultValue        = "0.05";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "fdthreshold";
    o->shorthand           = 'k';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->fdthreshold;
    o->cli_description     = "(optional) framewise displacement threshold";
    o->cli_arg_description = "<float>";
    o->defaultValue        = "1.0";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "mask";
    o->shorthand           = 'm';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->maskpath;
    o->cli_description     = "a mask specifying the region considered for the calculation of the DVARS. specify a brain/gm mask";
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
    o->name                = "modal";
    o->shorthand           = 'l';
    o->storage             = &p->useModal;
    o->cli_description     = "use the modal value for normalising the intensities as opposed to min/max";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
        
    // initialize the more advanced and rather unusual options
    o = rsUINewOption();
    o->name                = "dvars";
    o->shorthand           = 'd';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->dvarspath;
    o->cli_description     = "(optional) file where the DVARS values will be saved to";
    o->cli_arg_description = "<*.txt>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "fd";
    o->shorthand           = 'f';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->fdpath;
    o->cli_description     = "(optional) file to which the framewise displacement values will be saved to";
    o->cli_arg_description = "<*.txt>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "flagged";
    o->shorthand           = 'e';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->flaggedpath;
    o->cli_description     = "(optional) file to which the indices of all flagged frames will be saved to";
    o->cli_arg_description = "<*.txt>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
}

gboolean rsMotionScrubbingParseRPFormat(const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
    rsMotionScrubbingParameters *p = (rsMotionScrubbingParameters*) data;

    // parse value
    if ( ! strcmp(value, "fsl") ) {
        p->rpformat = RSTOOLS_REALIGNMENT_PARAMETER_FORMAT_FSL;
    } else if ( ! strcmp(value, "spm") ) {
        p->rpformat = RSTOOLS_REALIGNMENT_PARAMETER_FORMAT_SPM;
    } else {    
        // any other value should lead to an error
        g_set_error(
            error,
            G_OPTION_ERROR,
            G_OPTION_ERROR_BAD_VALUE,
            "%s: %s",
            option_name,
            "accepted values are 'fsl' and 'spm'"
        );
        return FALSE;
    }
    
    return TRUE;
}
