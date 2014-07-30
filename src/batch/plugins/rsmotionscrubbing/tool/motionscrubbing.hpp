#ifndef rstools_rsbatch_plugin_rsmotionscrubbing_tool_motionscrubbing_h
#define rstools_rsbatch_plugin_rsmotionscrubbing_tool_motionscrubbing_h

#include <iostream>
#include "src/batch/util/rstool.hpp"
#include "src/rsmotionscrubbing_common.h"
#include "src/rsmotionscrubbing_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsmotionscrubbing {
namespace tool {
    
class MotionScrubbing : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsMotionScrubbingParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rsmotionscrubbing::tool

#endif
