#ifndef rstools_rsbatch_plugin_rssmoothing_tool_bandpass_h
#define rstools_rsbatch_plugin_rssmoothing_tool_bandpass_h

#include <iostream>
#include "batch/util/rstool.hpp"
#include "rssmoothing_common.h"
#include "rssmoothing_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rssmoothing {
namespace tool {
    
class Smoothing : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsSmoothingParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rssmoothing::tool

#endif
