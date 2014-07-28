#ifndef rstools_rsbatch_plugin_rsbandpass_tool_bandpass_h
#define rstools_rsbatch_plugin_rsbandpass_tool_bandpass_h

#include <iostream>
#include "src/batch/util/rstool.hpp"
#include "src/rsregression_common.h"
#include "src/rsregression_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsregression {
namespace tool {
    
class Regression : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsRegressionParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rsregression::tool

#endif
