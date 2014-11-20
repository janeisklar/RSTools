#ifndef rstools_resample_common_h
#define rstools_resample_common_h

#include "rsresample_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsResampleInit(rsResampleParameters *p);
void rsResampleRun(rsResampleParameters *p);
void rsResampleDestroy(rsResampleParameters *p);
void rsResampleLanczosConvolve(double* signalOut, const double* signalIn, const int nVolsIn, const int nVolsOut, const int order, const double scaling);

#ifdef __cplusplus
}
#endif

#endif
