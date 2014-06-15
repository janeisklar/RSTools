#ifndef rstools_rsregression_common_h
#define rstools_rsregression_common_h

#include "rsregression_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsRegressionInit(rsRegressionParameters* p);
void rsRegressionRun(rsRegressionParameters *p);
void rsRegressionDestroy(rsRegressionParameters* p);

#ifdef __cplusplus
}
#endif

#endif
