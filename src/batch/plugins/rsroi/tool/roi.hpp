#ifndef rstools_rsbatch_plugin_rsroi_tool_roi_h
#define rstools_rsbatch_plugin_rsroi_tool_roi_h

#include <iostream>
#include "batch/util/rstool.hpp"
#include "rsroi_common.h"
#include "rsroi_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsroi {
namespace tool {
    
class Roi : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsRoiParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rsroi::tool

#endif
