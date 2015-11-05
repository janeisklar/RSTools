#ifndef rstools_deoblique_common_h
#define rstools_deoblique_common_h

#include "rsdeoblique_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsDeobliqueInit(rsDeobliqueParameters *p);
void rsDeobliqueRun(rsDeobliqueParameters *p);
void rsDeobliqueDestroy(rsDeobliqueParameters *p);

#ifdef __cplusplus
}
#endif

#endif
