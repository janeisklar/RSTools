#ifndef rstools_rsbatch_plugin_unix_tool_unix_h
#define rstools_rsbatch_plugin_unix_tool_unix_h

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "batch/util/rsunixtool.hpp"
#include "../task/unix.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {
namespace tool {
    
class Unix : public RSUnixTool {

public:
    void destroy();
    rsUIInterface* createUI();
    void printCallString(FILE *stream);
    
protected:
    void _init();
};

}}}}} // namespace rstools::batch::plugins::unix::tool

#endif
