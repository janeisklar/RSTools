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
    p->verbose              = FALSE;
    p->input                = NULL;
    p->output               = NULL;
    p->parametersValid      = FALSE;
    p->mask                 = NULL;
    p->context              = NULL;
    p->rp                   = NULL;
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
    g_option_context_free(p->context);
    rsFree(p);
}

void rsMotionScrubbingPrintHelp(rsMotionScrubbingParameters *p) {
    fprintf(stdout, "%s\n", g_option_context_get_help(p->context, TRUE, NULL));
}

rsMotionScrubbingParameters *rsMotionScrubbingParseParams(int argc, char * argv[]) {

    rsMotionScrubbingParameters *p = rsMotionScrubbingInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    // initialize the most common options
    GError *error = NULL;
    p->context = g_option_context_new("\n\nGiven a 4D-Nifti and a txt-file containing the 6 head realignment parameters, this application performs motion-scrubbing as described in: Power, Jonathan D., et al. \"Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion.\" Neuroimage 59.3 (2012): 2142-2154. APA");
    g_option_context_set_summary(p->context, RSTOOLS_VERSION_LABEL);

    /* long, short, flags, arg, arg_data, desc, arg_desc */
    GOptionEntry entries[] = {
      { "input",           'i', 0, G_OPTION_ARG_FILENAME, &p->inputpath,       "the 4D volume to be scrubbed", "<volume>" },
      { "output",          'o', 0, G_OPTION_ARG_FILENAME, &p->outputpath,      "the volume in which the result will be saved", "<volume>" },
      { "rp",              'r', 0, G_OPTION_ARG_FILENAME, &p->realignmentpath, "the file containing the realignment parameters", "<*.txt>" },
      { "dvars",           'd', 0, G_OPTION_ARG_FILENAME, &p->dvarspath,       "(optional) file where the DVARS values will be saved to", "<*.txt>" },
      { "fd",              'f', 0, G_OPTION_ARG_FILENAME, &p->fdpath,          "(optional) file to which the framewise displacement will be saved to", "<*.txt>" },
      { "flagged",         'e', 0, G_OPTION_ARG_FILENAME, &p->flaggedpath,     "(optional) file to which the indices of all flagged frames will be saved to", "<*.txt>" },
      { "dvarsthreshold",  'a', 0, G_OPTION_ARG_DOUBLE,   &p->dvarsthreshold,  "(optional) DVARs threshold", "<float>" },
      { "fdthreshold",     'k', 0, G_OPTION_ARG_DOUBLE,   &p->fdthreshold,     "(optional) framewise displacement threshold", "<float>" },
      { "mask",            'm', 0, G_OPTION_ARG_FILENAME, &p->maskpath,        "a mask specifying the region that is used for the computation of DVARs", "<volume>" },
      { "threads",         't', 0, G_OPTION_ARG_INT,      &p->threads,         "number of threads used for processing", "<n>" },
      { "verbose",         'v', 0, G_OPTION_ARG_NONE,     &p->verbose,         "show debug information", NULL},
      { NULL }
    };
    
    g_option_context_add_main_entries(p->context, entries, GETTEXT_PACKAGE);
    
    // check if parameters are valid
    if ( ! g_option_context_parse(p->context, &argc, &argv, &error) ) {
        fprintf(stderr, "option parsing failed: %s\n", error->message);
        return p;
    }
    
    p->parametersValid = TRUE;
    return p;
}
