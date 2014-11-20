#include "rsresample_ui.h"

rsResampleParameters *rsResampleInitParameters()
{
    rsResampleParameters *p = (rsResampleParameters*)rsMalloc(sizeof(rsResampleParameters));
    
    char *inputpath;
    char *outputpath;
    char *callString;
    
    double inputTR;
    double outputTR;
    
    BOOL verbose;
    BOOL parametersValid;
    
    rsNiftiFile *input;
    rsNiftiFile *output;

    rsUIInterface *interface;

    int threads;
    
    rsReportProgressCallback *progressCallback;
    
    
    p->inputpath        = NULL;
    p->outputpath       = NULL;
    p->callString       = NULL;
    p->inputTR          = -1.0;
    p->outputTR         = -1.0;
    p->verbose          = FALSE;
    p->input            = NULL;
    p->output           = NULL;
    p->order            = 2;
    p->regressors       = NULL;
    p->regressorInputs  = NULL;
    p->regressorOutputs = NULL;
    p->nRegressors      = 0;
    p->parametersValid  = FALSE;
    p->threads          = 1;
    p->interface        = NULL;
    p->progressCallback = NULL;
    
    return p;
}

void rsResampleFreeParams(rsResampleParameters *p)
{
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->input);
    rsFree(p->output);
    rsFree(p->callString);

    for (int i=0; i<p->nRegressors; i++) {
        rsFree(p->regressors[i]);
        rsFree(p->regressorInputs[i]);
        rsFree(p->regressorOutputs[i]);
    }
    
    rsFree(p->regressorInputs);
    rsFree(p->regressorOutputs);
    rsFree(p->regressors);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsResampleParameters *rsResampleParseParams(int argc, char * argv[])
{

    rsResampleParameters *p = rsResampleInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsResampleBuildInterface(p);
    
    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified(--input)!\n");
        return p;
    }
    
    if ( p->outputpath == NULL ) {
        fprintf(stderr, "An output path for the resampled data must be specified(--output)!\n");
        return p;
    }
    
    if ( p->inputTR < 0 || p->outputTR < 0 ) {
        fprintf(stderr, "Input and output sampling rates have to be specified (--TRin, --TRout)!\n");
        return p;
    }
    
    // parse additional regressor files that are supplied using -r or --regressor
    if ( p->regressors != NULL ) {
        p->nRegressors = g_strv_length(p->regressors);
        p->regressorInputs = (char**)rsMalloc(sizeof(char*)*((size_t)p->nRegressors));
        p->regressorOutputs = (char**)rsMalloc(sizeof(char*)*((size_t)p->nRegressors));
    
        // init with NULL so it isn't free'd in case of an error
        for(int i = 0; i < p->nRegressors; i++){
            p->regressorInputs[i] = NULL;
            p->regressorOutputs[i] = NULL;
        }
    
        // parse
        for(int i = 0; i < p->nRegressors; i++){
            const char *regressor = p->regressors[i];
            char *regressorIn;
            char *regressorOut;
                    
            if ( ! rsResampleParseRegressorPath(&regressorIn, &regressorOut, p->regressors[i]) ) {
                fprintf(stderr, "failed parsing regressor path: %s\n", p->regressors[i]);
                return p;
            }

            p->regressorInputs[i]  = regressorIn;
            p->regressorOutputs[i] = regressorOut;
        }
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsResampleBuildInterface(rsResampleParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Resamples a 4D-Nifti in the temporal domain using a lanczos filter.";
    
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the input volume that is going to be resampled";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "output";
    o->shorthand           = 'o';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->outputpath;
    o->cli_description     = "the output volume which is going to be resampled";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "TRin";
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->inputTR;
    o->cli_description     = "the sampling rate of the input data";
    o->cli_arg_description = "<TR in s>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "TRout";
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->outputTR;
    o->cli_description     = "the sampling rate of the output data";
    o->cli_arg_description = "<TR in s>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "regressor";
    o->shorthand           = 'r';
    o->type                = G_OPTION_ARG_STRING_ARRAY;
    o->storage             = &p->regressors;
    o->cli_description     = "can be used multiple times to specify a tabbed regressor file (may contain several columns). The input and output file should be specified together and separated by comma.";
    o->cli_arg_description = "<input,output>";
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
    o->name                = "order";
    o->shorthand           = 'a';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->order;
    o->cli_description     = "Order of the lanczos filter";
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
}

BOOL rsResampleParseRegressorPath(char **input, char **output, const char *arg)
{
    char *argCopy = (char*)rsMalloc((1+strlen(arg))*sizeof(char));
    sprintf(argCopy, "%s", arg);
    
    char *in  = strtok(argCopy, ",");
    char *out = strtok(NULL, ",");
        
    if ( strtok(NULL,",") != NULL || in == NULL || out == NULL ) {

        rsFree(argCopy);
        return FALSE;
    }
    
    *input  = rsString(in);
    *output = rsString(out);
    
    return TRUE;
}
