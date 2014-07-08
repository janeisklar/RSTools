#ifndef rstools_rsbatch_execution_bandpass_h
#define rstools_rsbatch_execution_bandpass_h

#include <iostream>
#include "tool.hpp"
#include "src/rsbandpass_common.h"

namespace rstools {
namespace batch {
namespace execution {
    
class Bandpass : public Tool {

public:
    void destroy();
    bool isEverythingFine();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsBandpassParameters *params;
};

}}} // namespace rstools::batch::execution

#endif
