#ifndef rstools_rsbatch_plugin_unix_tool_unix_h
#define rstools_rsbatch_plugin_unix_tool_unix_h

#include <iostream>
#include "src/batch/util/rstool.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {
namespace tool {
    
class Unix : public RSTool {

public:
    void destroy();
    bool isEverythingFine();
    rsUIInterface* createUI();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
};

}}}}} // namespace rstools::batch::plugins::unix::tool

#endif
