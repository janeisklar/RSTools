#ifndef rstools_rsbatch_plugin_rsscrubbing_tool_scrubbing_h
#define rstools_rsbatch_plugin_rsscrubbing_tool_scrubbing_h

#include <iostream>
#include "batch/util/rstool.hpp"
#include "rsscrubbing_common.h"
#include "rsscrubbing_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsscrubbing {
namespace tool {
    
class Scrubbing : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsScrubbingParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rsscrubbing::tool

#endif
