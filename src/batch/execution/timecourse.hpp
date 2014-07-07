#ifndef rstools_rsbatch_execution_timecourse_h
#define rstools_rsbatch_execution_timecourse_h

#include <iostream>
#include "tool.hpp"
#include "src/rstimecourse_common.h"

namespace rstools {
namespace batch {
namespace execution {
	
class Timecourse : public Tool {

public:
	void destroy();
	bool isEverythingFine();
	
protected:
	void _parseParams(int argc, char * argv[]);
	void _init();
	void _run();
	rsTimecourseParameters *params;
};

}}} // namespace rstools::batch::execution

#endif
