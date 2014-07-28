#ifndef rstools_rsbatch_plugin_rscorrelation_tool_correlation_h
#define rstools_rsbatch_plugin_rscorrelation_tool_correlation_h

#include <iostream>
#include "src/batch/util/rstool.hpp"
#include "src/rscorrelation_common.h"
#include "src/rscorrelation_ui.h"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rscorrelation {
namespace tool {
    
class Correlation : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsCorrelationParameters *params;
};

}}}}} // namespace rstools::batch::plugins::rscorrelation::tool

#endif
