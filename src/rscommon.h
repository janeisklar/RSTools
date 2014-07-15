#include "config.h"

#ifndef rstools_rscommon_h
#define rstools_rscommon_h

#include <stdio.h>
#include <nifti1.h>
#include <fslio.h>
#include <glib.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"

#define RSTOOLS_XSTR(x) #x
#define RSTOOLS_STR(x) RSTOOLS_XSTR(x)
#define RSTOOLS_VERSION_HASH_STRING RSTOOLS_VERSION_HASH
#define RSTOOLS_VERSION_DATE_STRING RSTOOLS_VERSION_DATE
#define RSTOOLS_VERSION VERSION " (" RSTOOLS_VERSION_HASH_STRING " - " RSTOOLS_VERSION_DATE_STRING ")"
#define RSTOOLS_VERSION_LABEL "RSTools " RSTOOLS_VERSION
#define RSTOOLS_DATA_DIR DATA_PATH
#define RSTOOLS_PLUGINS_DIR RSTOOLS_DATA_DIR "/" PACKAGE "/plugins"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double percentage;
    int run;
} rsReportProgressEvent;

typedef void (*rsReportProgressCallback_t)(rsReportProgressEvent *event, void *userdata);

typedef struct {
    rsReportProgressCallback_t cb;
    void *data;
} rsReportProgressCallback;

#ifdef __cplusplus
}
#endif

#endif
