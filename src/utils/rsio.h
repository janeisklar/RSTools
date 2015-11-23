#ifndef rstools_utils_io_h
#define rstools_utils_io_h

#include "rscommon.h"

// rare string to denote the end of a list of files
#define RSIO_LASTFILE "%%LAST--FILE%%"

#ifdef __cplusplus
extern "C" {
#endif

BOOL rsFileIsReadable(const char *path);
BOOL rsFileIsWritable(const char *path);
BOOL rsCheckInputs(const char **paths);
BOOL rsCheckOutputs(const char **paths);
BOOL rsReadline(FILE *f, char *line, int *length);

#ifdef __cplusplus
}
#endif

#endif
