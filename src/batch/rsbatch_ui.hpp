#ifndef rstools_rsbatch_ui_h
#define rstools_rsbatch_ui_h

#include <stdio.h>
#include <strings.h>
#include <glib.h>
#include "nifti/rsniftiutils.h"
#include "batch/util/rsjobparser.hpp"

#ifdef __cplusplus
extern "C" {
#endif

using namespace rstools::batch::util;

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
    int *tasksToSkip;
} rsBatchParameters;

rsBatchParameters* rsBatchParseParams(int argc, char * argv[]);
rsBatchParameters* rsBatchInitParameters();
void rsBatchFreeParams(rsBatchParameters* p);
void rsBatchPrintHelp(rsBatchParameters* p);
rsArgument* rsBatchParseArgument(char *arg);
gboolean rsBatchParseSkipParam(const gchar *option_name, const gchar *value, gpointer data, GError **error);

#ifdef __cplusplus
}
#endif

#endif
