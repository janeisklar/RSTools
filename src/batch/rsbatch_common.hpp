#ifndef rstools_rsbatch_common_h
#define rstools_rsbatch_common_h

#include "rsbatch_ui.hpp"
#include "execution/tool.hpp"
#include <sstream>
#include "src/utils/rsstring.h"

using namespace rstools::batch;

void rsBatchInit(rsBatchParameters* p);
void rsBatchRun(rsBatchParameters *p);
void rsBatchDestroy(rsBatchParameters* p);
void rsBatchPrintParameter(rsBatchParameters *p);
void rsBatchPrintExecutionError(execution::Tool *tool, int taskNum, char const * state);
void rsBatchShowJobOverview(rsBatchParameters* p, execution::Tool** tools);

#endif
