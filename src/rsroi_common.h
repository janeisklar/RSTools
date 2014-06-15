#ifndef rstools_rsroi_common_h
#define rstools_rsroi_common_h

#include "rsroi_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsRoiInit(rsRoiParameters *p);
void rsRoiRun(rsRoiParameters *p);
void rsRoiDestroy(rsRoiParameters *p);

#ifdef __cplusplus
}
#endif

#endif
