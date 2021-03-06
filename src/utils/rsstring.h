#ifndef rstools_utils_string_h
#define rstools_utils_string_h

#include "rscommon.h"
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "rsmemory.h"

#ifdef __cplusplus
extern "C" {
#endif

void  rsStringWordWrap(const char* inputString, char*** lineArray, size_t* nLines, const unsigned int lineLength);
char* rsStringConcat(const char *first, ...);
int   rsStringCompareCaseInsensitive(char const *a, char const *b);
char* rsTrimString(char *s);
char* rsLeftTrimString(char *s);
char* rsRightTrimString(char *s);
char* rsString(const char *s);
BOOL rsStringStartsWith(const char *str, const char *start);
BOOL rsStringEndsWith(const char *str, const char *end);
void rsStringAppend(char *str, const char *suffix);

#ifdef __cplusplus
}
#endif

#endif
