#ifndef rstools_zeropadding_common_h
#define rstools_zeropadding_common_h

#include "rszeropadding_ui.h"

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
} rsZeropaddingTransformations;

void rsZeropaddingInit(rsZeropaddingParameters *p);
void rsZeropaddingRun(rsZeropaddingParameters *p);
void rsZeropaddingDestroy(rsZeropaddingParameters *p);

#ifdef __cplusplus
}
#endif

#endif
