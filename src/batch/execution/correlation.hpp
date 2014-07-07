#ifndef rstools_rsbatch_execution_correlation_h
#define rstools_rsbatch_execution_correlation_h

#include <iostream>
#include "tool.hpp"
#include "src/rscorrelation_common.h"

namespace rstools {
namespace batch {
namespace execution {
	
class Correlation : public Tool {

public:
	void destroy();
	bool isEverythingFine();
	
protected:
	void _parseParams(int argc, char * argv[]);
	void _init();
	void _run();
	rsCorrelationParameters *params;
};

}}} // namespace rstools::batch::execution

#endif
