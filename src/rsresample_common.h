#ifndef rstools_resample_common_h
#define rstools_resample_common_h

#include "rsresample_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsResampleInit(rsResampleParameters *p);
void rsResampleRun(rsResampleParameters *p);
void rsResampleDestroy(rsResampleParameters *p);

#ifdef __cplusplus
}
#endif

#endif
