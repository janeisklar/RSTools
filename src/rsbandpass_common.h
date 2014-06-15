#ifndef rstools_rsbandpass_common_h
#define rstools_rsbandpass_common_h

#include "rsbandpass_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsBandpassInit(rsBandpassParameters *p);
void rsBandpassRun(rsBandpassParameters *p);
void rsBandpassDestroy(rsBandpassParameters *p);

#ifdef __cplusplus
}
#endif

#endif
