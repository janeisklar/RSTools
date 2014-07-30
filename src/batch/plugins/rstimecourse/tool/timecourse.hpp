#ifndef rstools_rsbatch_plugin_rstimecourse_tool_timecourse_h
#define rstools_rsbatch_plugin_rstimecourse_tool_timecourse_h

#include <iostream>
#include "src/batch/util/rstool.hpp"
#include "src/rstimecourse_common.h"
#include "src/rstimecourse_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rstimecourse {
namespace tool {
    
class Timecourse : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsTimecourseParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rstimecourse::tool

#endif
