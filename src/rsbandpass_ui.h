#ifndef rstools_bandpass_ui_h
#define  rstools_bandpass_ui_h

#include <stdio.h>
#include "nifti/rsniftiutils.h"
#include "maths/rsmathutils.h"
#include "utils/rsui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *inputpath;
    char *maskpath;
    char *savemaskpath;
    char *saveFilteredPath;
    char *saveAttenuationPath;
    char *callString;
    
    double freqLow;
    double freqHigh;
    double TR;
    
    BOOL verbose;
    BOOL parametersValid;
    
    BOOL keepMean;
    BOOL fftw;
    
    rsNiftiFile *input;
    rsNiftiFile *filteredOutput;
    double ***mask;

    rsUIInterface *interface;

    int threads;
    size_t wordsize;
    
    rsFFTFilterParams *fftParams;
    double rolloff;
    int rolloff_method;
    long paddedT;
    
    rsReportProgressCallback *progressCallback;
    
} rsBandpassParameters;

rsBandpassParameters *rsBandpassParseParams(int argc, char * argv[]);
rsBandpassParameters *rsBandpassInitParameters();
void rsBandpassFreeParams(rsBandpassParameters *p);
void rsBandpassBuildInterface(rsBandpassParameters *p);
    
#ifdef __cplusplus
}
#endif

#endif
