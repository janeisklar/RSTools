#ifndef rstools_rsbatch_execution_motionscrubbing_h
#define rstools_rsbatch_execution_motionscrubbing_h

#include <iostream>
#include "tool.hpp"
#include "src/rsmotionscrubbing_common.h"

namespace rstools {
namespace batch {
namespace execution {
    
class MotionScrubbing : public Tool {

public:
    void destroy();
    bool isEverythingFine();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsMotionScrubbingParameters *params;
};

}}} // namespace rstools::batch::execution

#endif
