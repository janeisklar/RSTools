#ifndef rstools_rsbatch_execution_unix_h
#define rstools_rsbatch_execution_unix_h

#include <iostream>
#include "tool.hpp"
#include <stdio.h>
#include <stdlib.h> 

namespace rstools {
namespace batch {
namespace execution {
    
class Unix : public Tool {

public:
    void destroy();
    bool isEverythingFine();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    bool executionSuccessful;
};

}}} // namespace rstools::batch::execution

#endif
