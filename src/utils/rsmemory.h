#ifndef rstools_utils_memory_h
#define rstools_utils_memory_h

#include "src/rscommon.h"

#define rsFree(x) { free(x); x=NULL; }

#ifdef __cplusplus
extern "C" {
#endif

void *rsMalloc(size_t size);

#ifdef __cplusplus
}
#endif

#endif
