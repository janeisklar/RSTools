#ifndef rstools_rsbatch_execution_regression_h
#define rstools_rsbatch_execution_regression_h

#include <iostream>
#include "tool.hpp"
#include "src/rsregression_common.h"

namespace rstools {
namespace batch {
namespace execution {
    
class Regression : public Tool {

public:
    void destroy();
    bool isEverythingFine();
    
protected:
    void _parseParams(int argc, char * argv[]);
    void _init();
    void _run();
    rsRegressionParameters *params;
};

}}} // namespace rstools::batch::execution

#endif
