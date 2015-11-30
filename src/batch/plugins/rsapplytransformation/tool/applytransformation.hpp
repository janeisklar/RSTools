#ifndef rstools_rsbatch_plugin_rsapplytransformation_tool_applytransformation_h
#define rstools_rsbatch_plugin_rsapplytransformation_tool_applytransformation_h

#include <iostream>
#include "batch/util/rstool.hpp"
#include "rsapplytransformation_common.h"
#include "rsapplytransformation_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsapplytransformation {
namespace tool {
    
class ApplyTransformation : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsApplyTransformationParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rsapplytransformation::tool

#endif
