#ifndef rstools_rsbatch_execution_roi_h
#define rstools_rsbatch_execution_roi_h

#include <iostream>
#include "tool.hpp"
#include "src/rsroi_common.h"

namespace rstools {
namespace batch {
namespace execution {
    
class Roi : public Tool {

public:
    void destroy();
    bool isEverythingFine();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsRoiParameters *params;
};

}}} // namespace rstools::batch::execution

#endif
