#ifndef rstools_rsbatch_plugin_rsorientation_tool_bandpass_h
#define rstools_rsbatch_plugin_rsorientation_tool_bandpass_h

#include <iostream>
#include "batch/util/rstool.hpp"
#include "rsorientation_common.h"
#include "rsorientation_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsorientation {
namespace tool {
    
class Orientation : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsOrientationParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rsorientation::tool

#endif
