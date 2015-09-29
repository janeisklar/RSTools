#ifndef rstools_rsinfo_common_h
#define rstools_rsinfo_common_h

#include "rsinfo_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsInfoInit(rsInfoParameters* p);
void rsInfoRun(rsInfoParameters *p);
void rsInfoDestroy(rsInfoParameters* p);

#ifdef __cplusplus
}
#endif

#endif
