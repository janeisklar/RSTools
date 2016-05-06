#ifndef rstools_rescale_common_h
#define rstools_rescale_common_h

#include "rsrescale_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsRescaleInit(rsRescaleParameters *p);
void rsRescaleRun(rsRescaleParameters *p);
void rsRescaleDestroy(rsRescaleParameters *p);

#ifdef __cplusplus
}
#endif

#endif
