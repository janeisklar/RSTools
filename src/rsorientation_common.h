#ifndef rstools_orientation_common_h
#define rstools_orientation_common_h

#include "rsorientation_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	BOOL xFlip;
	BOOL yFlip;
	BOOL zFlip;
	short xDim;
	short yDim;
	short zDim;
} rsOrientationTransformations;

void rsOrientationInit(rsOrientationParameters *p);
void rsOrientationRun(rsOrientationParameters *p);
void rsOrientationDestroy(rsOrientationParameters *p);

#ifdef __cplusplus
}
#endif

#endif
