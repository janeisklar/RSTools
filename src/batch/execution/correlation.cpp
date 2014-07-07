#include "correlation.hpp"

namespace rstools {
namespace batch {
namespace execution {
	
void Correlation::_parseParams(int argc, char * argv[])
{
	params = rsCorrelationParseParams(argc, argv);
}
	
void Correlation::_init()
{
	params->threads = this->threads;
	rsCorrelationInit(params);
}

void Correlation::_run()
{
	params->progressCallback = (rsReportProgressCallback*)rsMalloc(sizeof(rsReportProgressCallback));
	params->progressCallback->cb = (rsReportProgressCallback_t) Tool::showProgressCallback;
	params->progressCallback->data = (void*)oc;
		
	rsCorrelationRun(params);
	
	rsFree(params->progressCallback);
}

void Correlation::destroy()
{
	rsCorrelationDestroy(params);
}

bool Correlation::isEverythingFine()
{
	return params != NULL && params->parametersValid;
}

}}} // namespace rstools::batch::execution
