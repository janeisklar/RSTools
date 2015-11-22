#ifndef rstools_applytransformation_common_h
#define rstools_applytransformation_common_h

#include "rsapplytransformation_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsApplyTransformationInit(rsApplyTransformationParameters *p);
void rsApplyTransformationRun(rsApplyTransformationParameters *p);
void rsApplyTransformationDestroy(rsApplyTransformationParameters *p);

#ifdef __cplusplus
}
#endif

#endif
