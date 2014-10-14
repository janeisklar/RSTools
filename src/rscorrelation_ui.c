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
    p->interface             = NULL;
    p->progressCallback      = NULL;

    return p;
}

rsCorrelationParameters* rsCorrelationParseParams(int argc, char * argv[])
{
    
    rsCorrelationParameters *p = rsCorrelationInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    rsCorrelationBuildInterface(p);
    
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
        fprintf(stderr, "No output volume specified!\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsCorrelationBuildInterface(rsCorrelationParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "This tool will correlate the timecourse of every voxel in a 4D-nifti volume with a given timecourse.";
    p->interface->helpIndent    = 33;
    
    GOptionArgFunc cbConversion = (GOptionArgFunc)rsCorrelationParseConversionMode;
    GOptionArgFunc cbMonteCarlo = rsCorrelationParseMonteCarloParams;
    
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the volume for which the correlation of the timecourse for each voxel is computed";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "output";
    o->shorthand           = 'o';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->outputpath;
    o->cli_description     = "the volume in which the correlation values will be saved in";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "regressor";
    o->shorthand           = 'r';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->regressorpath;
    o->cli_description     = "txt-file containg a reference timecourse for which the correlation to the rest of the brain will be computed. if ommitted, the timecourse is expected to be supplied through the standard input instead";
    o->gui_description     = "txt-file containg a reference timecourse for which the correlation to the rest of the brain will be computed.";
    o->cli_arg_description = "<txt-file>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "conversion";
    o->shorthand           = 'c';
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbConversion;
    o->cli_description     = "<mode> specifies what is stored in the output file. It can take the follwing values:";
    o->cli_arg_description = "<mode>";
    o->defaultValue        = "z";
    rsUIOptionValue allowedValues[] = {
      {"none", "the correlation coefficients will be stored without converting them"},
      {"z",    "the correlation coefficients will be converted to z-values before being stored(default)"},
      {"t",    "the correlation coefficients will be converted to values of the T-statistic"},
      NULL
    };
    rsUISetOptionValues(o, allowedValues);
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "mask";
    o->shorthand           = 'm';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->maskpath;
    o->cli_description     = "a mask specifying the region that the correlation is perforned on (may be specified for improved performance)";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "comment";
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->commentpath;
    o->cli_description     = "adds a comment about the origin of thereference timecourse to the nifti header of the correlation map";
    o->cli_arg_description = "<txt>";    
    o->showInGUI           = FALSE;
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
        
    // initialize the more advanced and rather unusual options
    o = rsUINewOption();
    o->name                = "montecarlo";
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbMonteCarlo;
    o->cli_description     = "repeats the computation of the correlation n times and uses m randomly drawn timepoints in each run. eventually the average is being saved. (using it enforces the conversion to z-scores)";
    o->cli_arg_description = "<n>,<m>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "delay";
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->delay;
    o->cli_description     = "delay the regressor by <n> volumes(<n> * TR)";
    o->cli_arg_description = "<n>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
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
    rsUIDestroyInterface(p->interface);
    rsFree(p);
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
