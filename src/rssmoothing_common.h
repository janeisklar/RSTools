#ifndef rstools_scrubbing_common_h
#define rstools_scrubbing_common_h

#include "rssmoothing_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsSmoothingInit(rsSmoothingParameters *p);
void rsSmoothingRun(rsSmoothingParameters *p);
void rsSmoothingDestroy(rsSmoothingParameters *p);

#ifdef __cplusplus
}
#endif

#endif
