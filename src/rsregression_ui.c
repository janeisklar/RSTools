#include "rsregression_ui.h"

rsRegressionParameters* rsRegressionInitParameters()
{
    rsRegressionParameters *p = (rsRegressionParameters*)rsMalloc(sizeof(rsRegressionParameters));
    
    p->inputpath            = NULL;
    p->maskpath             = NULL;
    p->regressorspath       = NULL;
    p->regressorCommentPath = NULL;
    p->regressorComment     = NULL;
    p->comment              = NULL;
    p->savemaskpath         = NULL;
    p->saveBetasPath        = NULL;
    p->saveResidualsPath    = NULL;
    p->saveFittedPath       = NULL;
    p->callString           = NULL;
    p->freqLow              = -1.0;
    p->freqHigh             = -1.0;
    p->TR                   = -1.0;
    p->verbose              = FALSE;
    p->filterActive         = FALSE;
    p->zScoreRegression     = FALSE;
    p->nRegressors          = 0;
    p->nRegressorValues     = 0;
    p->regressors           = NULL;
    p->parametersValid      = FALSE;
    p->input                = NULL;
    p->betas                = NULL;
    p->residuals            = NULL;
    p->fitted               = NULL;
    p->mask                 = NULL;
    p->nyquist_frequency    = 0.0;
    p->bin_width            = 0.0;
    p->nFrequencyBinsLow    = 0;
    p->nFrequencyBinsHigh   = 0;
    p->nFrequencyBins       = 0;
    p->nFrequencyRegressors = 0;
    p->frequencyBins        = 0;
    p->nAllRegressors       = 0;
    p->allRegressors        = 0;
    p->threads              = 1;
    p->progressCallback     = NULL;
    p->interface            = NULL;

    return p;
}

rsRegressionParameters* rsRegressionParseParams(int argc, char * argv[])
{

    rsRegressionParameters *p = rsRegressionInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsRegressionBuildInterface(p);
    
    // parse
    p->parametersValid = rsUIParse(p->interface, argc, argv, (void*)p);
    return p;
}

void rsRegressionBuildInterface(rsRegressionParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Given a 4D-Nifti and a txt file with regressors(in the columns), this tool will perform a multiple linear regression on it.";
    
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the input volume to be regressed";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "residuals";
    o->shorthand           = 'r';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->saveResidualsPath;
    o->cli_description     = "the volume in which the residuals will be saved";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "fitted";
    o->shorthand           = 'f';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->saveFittedPath;
    o->cli_description     = "the volume in which the fitted values will be saved";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "betas";
    o->shorthand           = 'b';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->saveBetasPath;
    o->cli_description     = "the volume in which the beta-coefficients will be saved";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
        
    o = rsUINewOption();
    o->name                = "regressors";
    o->shorthand           = 'x';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->regressorspath;
    o->cli_description     = "a tabbed/spaced textfile containing the regressors with the different regressors in the columns and time course in the rows. Decimal numbers may be formatted like this: 1.23e+45";
    o->cli_arg_description = "<txt>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "mask";
    o->shorthand           = 'm';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->maskpath;
    o->cli_description     = "a mask specifying the region that the regression is perforned on (may be specified for improved performance)";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "regressorComment";
    o->shorthand           = 'c';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->regressorCommentPath;
    o->cli_description     = "adds a comment about the origin of the reference timecourse to the nifti header of the correlation map";
    o->cli_arg_description = "<txt>";    
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "zscore";
    o->shorthand           = 'z';
    o->storage             = &p->zScoreRegression;
    o->cli_description     = "performs a z-score based linear regression, i.e. the mean and standard deviation are removed from the input data and the regressors.";
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
    o->name                = "savemask";
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->savemaskpath;
    o->cli_description     = "optional path where the rescaled mask that was specified with the mask option will be saved. The saved file with have the same dimensions as the input volume.";
    o->cli_arg_description = "<volume>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "f1";
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->freqLow;
    o->cli_description     = "(optional) the lower frequency of the additional bandpass filter";
    o->cli_arg_description = "<frequency in Hz>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "f2";
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->freqHigh;
    o->cli_description     = "(optional) the upper frequency of the additional bandpass filter";
    o->cli_arg_description = "<frequency in Hz>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "TR";
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->TR;
    o->cli_description     = "(optional) the time to repeat (1/sampling frequency) that is used for the bandpass filter";
    o->cli_arg_description = "<rate in s>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
}

void rsRegressionFreeParams(rsRegressionParameters* p)
{
    rsFree(p->inputpath);
    rsFree(p->maskpath);
    rsFree(p->savemaskpath);
    rsFree(p->saveBetasPath);
    rsFree(p->saveResidualsPath);
    rsFree(p->saveFittedPath);
    rsFree(p->regressorCommentPath);
    rsFree(p->regressorComment);
    rsFree(p->comment);
    rsFree(p->regressorspath);
    if ( p->regressors != NULL ) {
        rsFree(p->regressors[0]);
        rsFree(p->regressors);
    }
    if ( p->mask != NULL ) {
        rsFree(p->mask[0][0]);
        rsFree(p->mask[0]);
        rsFree(p->mask);
    }
    rsFree(p->input);
    rsFree(p->betas);
    rsFree(p->fitted);
    rsFree(p->frequencyBins);
    if ( p->filterActive ) {
        rsFree(p->allRegressors[0]);
        rsFree(p->allRegressors);
    }
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}
