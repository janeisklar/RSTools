#include "rsbandpass_ui.h"

rsBandpassParameters *rsBandpassInitParameters() {
    rsBandpassParameters *p = (rsBandpassParameters*)rsMalloc(sizeof(rsBandpassParameters));
    
    p->inputpath            = NULL;
    p->maskpath             = NULL;
    p->savemaskpath         = NULL;
    p->saveFilteredPath     = NULL;
    p->saveAttenuationPath  = NULL;
    p->callString           = NULL;
    p->fftParams            = NULL;
    p->paddedT              = 0;
    p->freqLow              = -1.0;
    p->freqHigh             = -1.0;
    p->TR                   = -1.0;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->filteredOutput       = NULL;
    p->parametersValid      = FALSE;
    p->mask                 = NULL;
    p->threads              = 1;
    p->rolloff_method       = RSFFTFILTER_CUTOFF;
    p->rolloff              = 10.0;
    p->keepMean             = FALSE;
    p->context              = NULL;
    p->progressCallback     = NULL;
    
    return p;
}

void rsBandpassFreeParams(rsBandpassParameters *p) {
    free(p->inputpath);
    free(p->maskpath);
    free(p->savemaskpath);
    free(p->saveFilteredPath);
    free(p->saveAttenuationPath);
    free(p->input);
    free(p->filteredOutput);
    free(p->callString);
    free(p->fftParams);
    g_option_context_free(p->context);
    free(p);
}

void rsBandpassPrintHelp(rsBandpassParameters *p) {
    fprintf(stdout, "%s\n", g_option_context_get_help(p->context, TRUE, NULL));
}

rsBandpassParameters *rsBandpassParseParams(int argc, char * argv[]) {

    rsBandpassParameters *p = rsBandpassInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    // initialize the most common options
    GError *error = NULL;
    p->context = g_option_context_new("\n\nGiven a 4D-Nifti and a frequency band this tool will apply FFT-based temporal filtering in the specified frequency range onto the dataset.");
    g_option_context_set_summary(p->context, RSTOOLS_VERSION_LABEL);

    /* long, short, flags, arg, arg_data, desc, arg_desc */
    GOptionEntry entries[] = {
      { "input",           'i', 0, G_OPTION_ARG_FILENAME, &p->inputpath,           "the input volume that is going to be filtered", "<volume>" },
      { "filtered",        'f', 0, G_OPTION_ARG_FILENAME, &p->saveFilteredPath,    "the output volume in which the filtered data will be saved", "<volume>" },
      { "f1",              'l', 0, G_OPTION_ARG_DOUBLE,   &p->freqLow,             "the lower frequency of the bandpass filter", "<frequency in Hz>" },
      { "f2",              'u', 0, G_OPTION_ARG_DOUBLE,   &p->freqHigh,            "the upper frequency of the bandpass filter", "<frequency in Hz>" },
      { "TR",              'r', 0, G_OPTION_ARG_DOUBLE,   &p->TR,                  "the time to repeat (1/sampling frequency) that is used for the FFT", "<rate in s>" },
      { "mask",            'm', 0, G_OPTION_ARG_FILENAME, &p->maskpath,            "a mask specifying the region that the filter is applied on (may be specified for improved performance)", "<volume>" },
      { "keepMean",        'k', 0, G_OPTION_ARG_NONE,     &p->keepMean,            "retains the first bin of the FFT (the mean) independent of it being included in the selected frequency range", NULL },
      { "threads",         't', 0, G_OPTION_ARG_INT,      &p->threads,             "number of threads used for processing", "<n>" },
      { "verbose",         'v', 0, G_OPTION_ARG_NONE,     &p->verbose,             "show debug information", NULL},
      { NULL }
    };
    
    g_option_context_add_main_entries(p->context, entries, GETTEXT_PACKAGE);
    
    // initialize the more advanced and rather unusual options
    GOptionGroup *g = g_option_group_new("extended", "Extended options", "Additional options that are rarely going to be used", NULL, NULL);
    g_option_context_add_group(p->context, g);
    
    GOptionEntry extended_entries[] = {
      { "savemask",          0, 0, G_OPTION_ARG_FILENAME, &p->savemaskpath,        "optional path where the rescaled mask specified with -mask will be saved. The saved file with have the same dimensions as the input volume.", "volume" },
#if RS_FFTW_ENABLED == 1
      { "fftw",              0, 0, G_OPTION_ARG_NONE,     &p->fftw,                "use FFTW3 instead of GSL for FFT", NULL },
#endif  
      { "saveattenuation",   0, 0, G_OPTION_ARG_FILENAME, &p->saveAttenuationPath, "save txt file that contains the bin's frequencies and the corresponding attenuation weight that was used.", "txt-file" },
      { "sigmoidrolloff",    0, 0, G_OPTION_ARG_DOUBLE,   &p->rolloff,             "uses a sigmoid function for rolling off the passband. The specified number controls how fast it is rolled off with higher numbers corresponding to a quicker rolloff. A good starting point would be 10, then double-check by saving the attenuation file. (does not work with FFTW3)", "double" },
      { NULL }
    };
    
    g_option_group_add_entries(g, extended_entries);
    
    // check if parameters are valid
    if ( ! g_option_context_parse(p->context, &argc, &argv, &error) ) {
        fprintf(stderr, "option parsing failed: %s\n", error->message);
        return p;
    }
    
    p->parametersValid = TRUE;
    return p;
}
