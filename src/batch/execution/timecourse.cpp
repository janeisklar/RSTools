#include "timecourse.hpp"

namespace rstools {
namespace batch {
namespace execution {
	
void Timecourse::_parseParams(int argc, char * argv[])
{
	params = rsTimecourseParseParams(argc, argv);
}
	
void Timecourse::_init()
{
	params->threads = this->threads;
	rsTimecourseInit(params);
}

void Timecourse::_run()
{
	rsTimecourseRun(params);
}

void Timecourse::destroy()
{
	rsTimecourseDestroy(params);
}

bool Timecourse::isEverythingFine()
{
	return params != NULL && params->parametersValid;
}

}}} // namespace rstools::batch::execution