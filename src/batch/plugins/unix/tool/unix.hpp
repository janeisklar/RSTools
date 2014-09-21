#ifndef rstools_rsbatch_plugin_unix_tool_unix_h
#define rstools_rsbatch_plugin_unix_tool_unix_h

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "batch/util/rstool.hpp"
#include "../task/unix.hpp"

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
    void printCallString(FILE *stream);
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    
    rstools::batch::plugins::unix::task::Unix* getUnixTask();
    
    bool executionSuccessful;
};

}}}}} // namespace rstools::batch::plugins::unix::tool

#endif
