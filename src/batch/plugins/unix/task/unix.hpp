#ifndef rstools_rsbatch_plugin_unix_task_unix_h
#define rstools_rsbatch_plugin_unix_task_unix_h

#include <iostream>
#include "batch/util/rsunixtask.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {
namespace task {
    
class Unix : public RSUnixTask {

public: 
    Unix(const char* code, const char* name);
    
    void fillInJobArguments(RSJob* job, RSJobParser* parser);    
    void parseTaskFromXml(DOMNodeIterator* walker, DOMNode* &current_node);    
    char* toXml();
    
    void setCmd(const char* cmd);
    char* getCmd();
    
protected:
};

}}}}} // namespace rstools::batch::plugins::unix::task

#endif
