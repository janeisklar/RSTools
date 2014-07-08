#include "rscorrelation_ui.h"

rsCorrelationParameters* rsCorrelationInitParameters()
{
    rsCorrelationParameters *p = (rsCorrelationParameters*)rsMalloc(sizeof(rsCorrelationParameters));
    
    p->inputpath             = NULL;
    p->regressorpath         = NULL;
    p->maskpath              = NULL;
    p->outputpath            = NULL;
    p->savemaskpath          = NULL;
    p->commentpath           = NULL;
    p->comment               = NULL;
    p->callString            = NULL;
    p->delay                 = 0;
    p->verbose               = FALSE;
    p->input                 = NULL;
    p->output                = NULL;
    p->correlation           = NULL;
    p->parametersValid       = FALSE;
    p->mask                  = NULL;
    p->nRegressorValues      = 0;
    p->regressor             = NULL;
    p->threads               = 1;
    p->conversionMode        = RSTOOLS_CORRELATION_CONVERSION_Z;
    p->monteCarloRepetitions = 0;
    p->monteCarloSampleSize  = 0;
    p->context               = NULL;
    p->progressCallback      = NULL;

    return p;
}

rsCorrelationParameters* rsCorrelationParseParams(int argc, char * argv[])
{
    
    rsCorrelationParameters *p = rsCorrelationInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    // initialize the most common options
    GError *error = NULL;
    p->context = g_option_context_new("\n\nThis tool will correlate the timecourse of of every voxel in a 4D-nifti volume with a timecourse that is supplied via standard input.");
    GOptionGroup *g = g_option_group_new("common", "Common options", "The most commonly used options", (void*)p, NULL);
    
    g_option_context_set_summary(p->context, RSTOOLS_VERSION_LABEL);
    
    GOptionArgFunc cbConversion = (GOptionArgFunc)rsCorrelationParseConversionMode;
    
     /* long, short, flags, arg, arg_data, desc, arg_desc */
    GOptionEntry entries[] = {
      { "input",           'i', 0, G_OPTION_ARG_FILENAME, &p->inputpath,     "the volume for which the correlation of the timecourse for each voxel is computed", "<volume>" },
      { "output",          'o', 0, G_OPTION_ARG_FILENAME, &p->outputpath,    "the volume in which the correlation values will be saved in", "<volume>" },
      { "regressor",       'r', 0, G_OPTION_ARG_FILENAME, &p->regressorpath, "txt-file containg a timecourse for which the correlation to the rest of the brain will be computed. if ommitted, the timecourse is expected to be supplied through the standard input instead", "<txt-file>" },
      { "conversion",      'c', 0, G_OPTION_ARG_CALLBACK, cbConversion,      "<mode> specifies what is stored in the output file, it can take the follwing values:\n\n\t'none' - the correlation coefficients will be stored without converting them\n\t'z'    - the correlation coefficients will be converted to z-values before being stored(default)\n\t't'    - the correlation coefficients will be converted to values of the T-statistic\n", "<mode>" },
      { "mask",            'm', 0, G_OPTION_ARG_FILENAME, &p->maskpath,      "a mask specifying the region that the correlation is perforned on (may be specified for improved performance)", "<volume>" },
      { "comment",           0, 0, G_OPTION_ARG_STRING,   &p->commentpath,   "adds a comment about the origin of thereference timecourse to the nifti header of the correlation map", "<txt-file>" },
      { "threads",         't', 0, G_OPTION_ARG_INT,      &p->threads,       "number of threads used for processing", "<n>" },
      { "verbose",         'v', 0, G_OPTION_ARG_NONE,     &p->verbose,       "show debug information", NULL},
      { NULL }
    };
    
    g_option_context_set_main_group(p->context, g);
    g_option_group_add_entries(g, entries);

    // initialize the more advanced and rather unusual options
    g = g_option_group_new("extended", "Extended options", "Additional options that are rarely going to be used", (void*)p, NULL);
    g_option_context_add_group(p->context, g);

    GOptionArgFunc cbMonteCarlo = rsCorrelationParseMonteCarloParams;

    GOptionEntry extended_entries[] = {
      { "montecarlo",        0, 0, G_OPTION_ARG_CALLBACK, &cbMonteCarlo,   "repeats the computation of the correlation n times and uses m randomly drawn timepoints in each run. eventually the average is being saved. (using it enforces the conversion to z-scores)", "<n>,<m>" },
      { "delay",           'd', 0, G_OPTION_ARG_INT,      &p->delay,       "delay the regressor by <n> volumes(<n> * TR)", "<n>" },
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


void rsCorrelationFreeParams(rsCorrelationParameters* p)
{
    rsFree(p->inputpath);
    rsFree(p->maskpath);
    rsFree(p->savemaskpath);
    rsFree(p->outputpath);
    rsFree(p->commentpath);
    rsFree(p->comment);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->correlation);
    rsFree(p->mask);
    rsFree(p->regressor);
    g_option_context_free(p->context);
    rsFree(p);
}

void rsCorrelationPrintHelp(rsCorrelationParameters* p)
{
    fprintf(stdout, "%s\n", g_option_context_get_help(p->context, TRUE, NULL));
}

gboolean rsCorrelationParseConversionMode(const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
    rsCorrelationParameters *p = (rsCorrelationParameters*) data;

    // parse value
    if ( ! strcmp(value, "none") ) {
        p->conversionMode = RSTOOLS_CORRELATION_CONVERSION_NONE;
        return TRUE;
    } else if ( ! strcmp(value, "z") ) {
        p->conversionMode = RSTOOLS_CORRELATION_CONVERSION_Z;
        return TRUE;
    } else if ( ! strcmp(value, "t") ) {
        p->conversionMode = RSTOOLS_CORRELATION_CONVERSION_T;
        return TRUE;
    } 
    
    // any other value should lead to an error
    g_set_error(
        error,
        G_OPTION_ERROR,
        G_OPTION_ERROR_BAD_VALUE,
        "%s: %s",
        option_name,
        "accepted values are 'none', 't' and 'z'"
    );
    
    return FALSE;
}

gboolean rsCorrelationParseMonteCarloParams(const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
    rsCorrelationParameters *p = (rsCorrelationParameters*) data;

    // copy value(const)
    size_t length = strlen(value);
    char v[length+1];
    sprintf(&v[0], "%s", value);

    // parse value
    char *repetitions;
    char *samples;
    repetitions = strtok(v, ",");
    samples     = strtok(NULL,  ",");

    p->monteCarloRepetitions = atoi(repetitions);
    p->monteCarloSampleSize  = atoi(samples);
    
    // if we were given exactly 2 numbers separated by a comma return success 
    if ( strtok(NULL,",") == NULL && p->monteCarloRepetitions > 0 && p->monteCarloSampleSize > 0 ) {
        return TRUE;
    }
    
    // anything else should lead to an error
    g_set_error(
        error,
        G_OPTION_ERROR,
        G_OPTION_ERROR_BAD_VALUE,
        "%s: %s",
        option_name,
        "format should be n,m where n is the number of repetitions and m the amount of samples that are drawn in each run"
    );
    
    return FALSE;
}
