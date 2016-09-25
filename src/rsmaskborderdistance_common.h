#ifndef rstools_rsmaskborderdistance_common_h
#define rstools_rsmaskborderdistance_common_h

#include "rsmaskborderdistance_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsMaskBorderDistanceInit(rsMaskBorderDistanceParameters *p);
void rsMaskBorderDistanceRun(rsMaskBorderDistanceParameters *p);
void rsMaskBorderDistanceDestroy(rsMaskBorderDistanceParameters *p);

#ifdef __cplusplus
}
#endif

#endif
