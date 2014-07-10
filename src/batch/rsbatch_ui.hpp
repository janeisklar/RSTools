#ifndef rstools_rsbatch_ui_h
#define rstools_rsbatch_ui_h

#include <stdio.h>
#include <strings.h>
#include <glib.h>
#include "src/nifti/rsniftiutils.h"
#include "rsjobparser.hpp"

typedef struct {
    char *jobpath;
    RSJob *job;
    rsArgument **arguments;
    short nArguments;
    BOOL verbose;
    BOOL quiet;
    BOOL showOverview;
    BOOL parametersValid;
    GOptionContext *context;
    char *viewArgument;
    int threads;
} rsBatchParameters;

rsBatchParameters* rsBatchParseParams(int argc, char * argv[]);
rsBatchParameters* rsBatchInitParameters();
void rsBatchFreeParams(rsBatchParameters* p);
void rsBatchPrintHelp(rsBatchParameters* p);
rsArgument* rsBatchParseArgument(char *arg);

#endif
