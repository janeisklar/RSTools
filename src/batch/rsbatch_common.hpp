#ifndef rstools_rsbatch_common_h
#define rstools_rsbatch_common_h

#include "rsbatch_ui.hpp"
#include "src/batch/util/rstool.hpp"
#include "src/batch/util/pluginmanager.hpp"
#include <sstream>
#include "src/utils/rsstring.h"

using namespace rstools::batch::util;

void rsBatchInit(rsBatchParameters* p);
void rsBatchRun(rsBatchParameters *p);
void rsBatchDestroy(rsBatchParameters* p);
void rsBatchPrintParameter(rsBatchParameters *p);
void rsBatchPrintExecutionError(RSTool *tool, int taskNum, char const * state);
void rsBatchShowJobOverview(rsBatchParameters* p, RSTool** tools);
BOOL rsBatchTaskShouldBeSkipped(rsBatchParameters* p, int taskId);

#endif
