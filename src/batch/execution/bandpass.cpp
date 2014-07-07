#include "bandpass.hpp"

namespace rstools {
namespace batch {
namespace execution {
	
void Bandpass::_parseParams(int argc, char * argv[])
{
	params = rsBandpassParseParams(argc, argv);
}
	
void Bandpass::_init()
{
	params->threads = this->threads;	
	rsBandpassInit(params);
}

void Bandpass::_run()
{
	params->progressCallback = (rsReportProgressCallback*)rsMalloc(sizeof(rsReportProgressCallback));
	params->progressCallback->cb = (rsReportProgressCallback_t) Tool::showProgressCallback;
	params->progressCallback->data = (void*)oc;
		
	rsBandpassRun(params);
	
	rsFree(params->progressCallback);
}

void Bandpass::destroy()
{
	rsBandpassDestroy(params);
}

bool Bandpass::isEverythingFine()
{
	return params != NULL && params->parametersValid;
}

}}} // namespace rstools::batch::execution
