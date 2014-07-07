#ifndef rstools_rsbatch_common_h
#define rstools_rsbatch_common_h

#include "rsbatch_ui.hpp"
#include "execution/timecourse.hpp"
#include "execution/regression.hpp"
#include "execution/bandpass.hpp"
#include "execution/motionscrubbing.hpp"
#include "execution/correlation.hpp"
#include "execution/roi.hpp"
#include "execution/unix.hpp"

using namespace rstools::batch;

void rsBatchInit(rsBatchParameters* p);
void rsBatchRun(rsBatchParameters *p);
void rsBatchDestroy(rsBatchParameters* p);
void rsBatchPrintParameter(rsBatchParameters *p);
void rsBatchPrintExecutionError(execution::Tool *tool, int taskNum, char const * state);

#endif
