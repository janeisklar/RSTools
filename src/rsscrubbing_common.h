#ifndef rstools_scrubbing_common_h
#define rstools_scrubbing_common_h

#include "rsscrubbing_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsScrubbingInit(rsScrubbingParameters *p);
void rsScrubbingRun(rsScrubbingParameters *p);
void rsScrubbingDestroy(rsScrubbingParameters *p);

BOOL rsLoadIndexVector(char *path, BOOL **indices, long length);

#ifdef __cplusplus
}
#endif

#endif
