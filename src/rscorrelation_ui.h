#ifndef rstools_rscorrelation_ui_h
#define rstools_rscorrelation_ui_h

#include <stdio.h>
#include <strings.h>
#include <nifti1.h>
#include <fslio.h>
#include <glib.h>
#include <dlfcn.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"

#define RSTOOLS_CORRELATION_CONVERSION_NONE 1
#define RSTOOLS_CORRELATION_CONVERSION_Z    2
#define RSTOOLS_CORRELATION_CONVERSION_T    3

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *inputpath;
    char *regressorpath;
    char *maskpath;
    char *outputpath;
    char *commentpath;
    char *comment;
    char *callString;
    char *savemaskpath;
    char *saveBetasPath;
    char *saveResidualsPath;
    char *saveFittedPath;
    
    short delay;
    
    int conversionMode;

    unsigned int monteCarloRepetitions;
    unsigned int monteCarloSampleSize;
    
    BOOL verbose;
    BOOL parametersValid;
    
    unsigned int nRegressorValues;
    double *regressor;
    
    rsNiftiFile *input;
    rsNiftiFile *output;
    double ***mask;
    double ***correlation;

    GOptionContext *context;
    
    int threads;
    size_t wordsize;

    rsReportProgressCallback *progressCallback;
    
} rsCorrelationParameters;

rsCorrelationParameters* rsCorrelationParseParams(int argc, char * argv[]);
rsCorrelationParameters* rsCorrelationInitParameters();
void rsCorrelationFreeParams(rsCorrelationParameters* p);
void rsCorrelationPrintHelp(rsCorrelationParameters* p);

gboolean rsCorrelationParseConversionMode(const gchar *option_name, const gchar *value, gpointer data, GError **error);
gboolean rsCorrelationParseMonteCarloParams(const gchar *option_name, const gchar *value, gpointer data, GError **error);

#ifdef __cplusplus
}
#endif

#endif
