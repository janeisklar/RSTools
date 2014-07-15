#ifndef rstools_utils_string_h
#define rstools_utils_string_h

#include "src/rscommon.h"
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "rsmemory.h"

#ifdef __cplusplus
extern "C" {
#endif

void  rsStringWordWrap(const char* inputString, char*** lineArray, size_t* nLines, const unsigned int lineLength);
char* rsStringConcat(char *first, ...);
int   rsStringCompareCaseInsensitive(char const *a, char const *b);

#ifdef __cplusplus
}
#endif

#endif
