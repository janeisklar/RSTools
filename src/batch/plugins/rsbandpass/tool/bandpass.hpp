#ifndef rstools_rsbatch_plugin_rsbandpass_tool_bandpass_h
#define rstools_rsbatch_plugin_rsbandpass_tool_bandpass_h

#include <iostream>
#include "batch/util/rstool.hpp"
#include "rsbandpass_common.h"
#include "rsbandpass_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsbandpass {
namespace tool {
    
class Bandpass : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsBandpassParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rsbandpass::tool

#endif
