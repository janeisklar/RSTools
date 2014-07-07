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
	p->context              = NULL;
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

    return p;
}

rsRegressionParameters* rsRegressionParseParams(int argc, char * argv[])
{
    rsRegressionParameters *p = rsRegressionInitParameters();
	p->callString = rsMergeStringArray(argc, argv);

	// initialize the most common options
	GError *error = NULL;
	p->context = g_option_context_new("\n\n Given a 4D-Nifti and a txt file with regressors(in the columns), this tool will perform a multiple linear regression on it.");
	GOptionGroup *g = g_option_group_new("common", "Common options", "The most commonly used options", (void*)p, NULL);
	
	g_option_context_set_summary(p->context, RSTOOLS_VERSION_LABEL);
	
	 /* long, short, flags, arg, arg_data, desc, arg_desc */
	GOptionEntry entries[] = {
	  { "input",            'i', 0, G_OPTION_ARG_FILENAME, &p->inputpath,            "the volume to be regressed", "<volume>" },
	  { "residuals",        'r', 0, G_OPTION_ARG_FILENAME, &p->saveResidualsPath,    "the volume in which the residuals will be saved", "<volume>" },
	  { "fitted",           'f', 0, G_OPTION_ARG_FILENAME, &p->saveFittedPath,       "the volume in which the fitted values will be saved", "<volume>" },
	  { "betas",            'b', 0, G_OPTION_ARG_FILENAME, &p->saveBetasPath,        "the volume in which the beta-coefficients will be saved", "<volume>" },
	  { "regressors",       'x', 0, G_OPTION_ARG_FILENAME, &p->regressorspath,       "a tabbed/spaced textfile containing the regressors with the different regressors in the columns and time course in the rows. Decimal numbers may be formatted like this: 1.23e+45", "<txt>" },
	  { "mask",             'm', 0, G_OPTION_ARG_FILENAME, &p->maskpath,             "a mask specifying the region that the regression is perforned on (may be specified for improved performance)", "<volume>" },
	  { "regressorComment", 'c', 0, G_OPTION_ARG_FILENAME, &p->regressorCommentPath, "adds a comment about the origin of thereference timecourse to the nifti header of the correlation map", "<txt>" },
	  { "zscore",           'z', 0, G_OPTION_ARG_NONE,     &p->zScoreRegression,     "performs a z-score based linear regression, i.e. the mean and standard deviation are removed from the input data and the regressors.", NULL },
	  { "threads",          't', 0, G_OPTION_ARG_INT,      &p->threads,              "number of threads used for processing", "<n>" },
	  { "verbose",          'v', 0, G_OPTION_ARG_NONE,     &p->verbose,              "show debug information", NULL},
	  { NULL }
	};
	
	g_option_context_set_main_group(p->context, g);
	g_option_group_add_entries(g, entries);

	// initialize the more advanced and rather unusual options
	g = g_option_group_new("extended", "Extended options", "Additional options that are rarely going to be used", (void*)p, NULL);
	g_option_context_add_group(p->context, g);

	GOptionEntry extended_entries[] = {
	  { "savemask",          0, 0, G_OPTION_ARG_FILENAME, &p->savemaskpath,        "optional path where the rescaled mask specified with -mask will be saved. The saved file with have the same dimensions as the input volume.", "volume" },
	  { "f1",                0, 0, G_OPTION_ARG_DOUBLE,   &p->freqLow,             "(optional) the lower frequency of the additional bandpass filter", "<frequency in Hz>" },
	  { "f2",                0, 0, G_OPTION_ARG_DOUBLE,   &p->freqHigh,            "(optional) the upper frequency of the additional bandpass filter", "<frequency in Hz>" },
	  { "TR",                0, 0, G_OPTION_ARG_DOUBLE,   &p->TR,                  "(optional) the time to repeat (1/sampling frequency) that is used for the bandpass filter", "<rate in s>" },
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
	g_option_context_free(p->context);
	rsFree(p);
}

void rsRegressionPrintHelp(rsRegressionParameters* p)
{
	fprintf(stdout, "%s\n", g_option_context_get_help(p->context, TRUE, NULL));
}
