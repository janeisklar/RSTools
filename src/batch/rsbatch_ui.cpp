#include "rsbatch_ui.hpp"

using namespace rstools::batch::util;

rsBatchParameters* rsBatchInitParameters()
{
    rsBatchParameters *p = (rsBatchParameters*)rsMalloc(sizeof(rsBatchParameters));
    
    p->jobpath               = NULL;
    p->arguments             = NULL;
    p->nArguments            = 0;
    p->context               = NULL;
    p->verbose               = FALSE;
    p->parametersValid       = FALSE;
    p->viewArgument          = NULL;
    p->showOverview          = FALSE;
    p->threads               = 1;
    p->quiet                 = FALSE;

    return p;
}

rsBatchParameters* rsBatchParseParams(int argc, char * argv[])
{    
    rsBatchParameters *p = rsBatchInitParameters();

    // initialize the most common options
    GError *error = NULL;
    p->context = g_option_context_new("\n\n...");
    GOptionGroup *g = g_option_group_new("common", "Common options", "The most commonly used options", (void*)p, NULL);
    
    g_option_context_set_summary(p->context, RSTOOLS_VERSION_LABEL);
    
    char **arguments = NULL;
    
     /* long, short, flags, arg, arg_data, desc, arg_desc */
    GOptionEntry entries[] = {
      { "job",           'j', 0, G_OPTION_ARG_FILENAME,     &p->jobpath,      "the jobfile that is to be executed", "<*.rsjob>" },
      { "argument",      'A', 0, G_OPTION_ARG_STRING_ARRAY, &arguments,       "additional arguments that will be passed onto the job-file, i.e. -A prefix=/usr/data/study2014/subjects", NULL },
      { "view",          'V', 0, G_OPTION_ARG_STRING,       &p->viewArgument, "prints the value for a the specified parameter in the job-file and quits the program. the job will not be executed", "argName"},
      { "overview",      'o', 0, G_OPTION_ARG_NONE,         &p->showOverview, "list tasks in the job and the corresponding parameters, but does not execute any task", NULL},
      { "threads",       't', 0, G_OPTION_ARG_INT,          &p->threads,      "number of threads used for processing", "<n>" },
      { "quiet",         'q', 0, G_OPTION_ARG_NONE,         &p->quiet,        "don't output the progress", NULL},
      { NULL }
    };
    
    g_option_context_set_main_group(p->context, g);
    g_option_group_add_entries(g, entries);

    // check if parameters are valid
    if ( ! g_option_context_parse(p->context, &argc, &argv, &error) ) {
        fprintf(stderr, "option parsing failed: %s\n", error->message);
        return p;
    }

    // parse the additional arguments that are supplied using -A
    if ( arguments != NULL ) {
        p->nArguments = g_strv_length(arguments);
        p->arguments = (rsArgument**)rsMalloc(sizeof(rsArgument*)*(size_t)p->nArguments);
    }
    
    for(int i = 0; i< p->nArguments; i++){
        char *argument = arguments[i];
        
        p->arguments[i] = rsBatchParseArgument(argument);
        
        if ( p->arguments[i] == NULL ) {
            fprintf(stderr, "failed parsing additional argument: %s\n", argument);
            return p;
        }       
    }

    p->parametersValid = TRUE;
    return p;
}


void rsBatchFreeParams(rsBatchParameters* p)
{
    free(p->jobpath);
    free(p->arguments);
    g_option_context_free(p->context);
    free(p);
}

void rsBatchPrintHelp(rsBatchParameters* p)
{
    fprintf(stdout, "%s\n", g_option_context_get_help(p->context, TRUE, NULL));
}

rsArgument* rsBatchParseArgument(char *arg)
{
    char *argCopy = (char*)rsMalloc((1+strlen(arg))*sizeof(char));
    sprintf(argCopy, "%s", arg);
    
    rsArgument *argument = (rsArgument*)rsMalloc(sizeof(rsArgument));
    argument->key   = strtok(argCopy, "=");
    argument->value = strtok(NULL, "=");
    
    if ( strtok(NULL,"=") != NULL || argument->key == NULL || argument->value == NULL ) {
        rsFree(argument);
        rsFree(argCopy);
        return NULL;
    }
    
    return argument;
}
