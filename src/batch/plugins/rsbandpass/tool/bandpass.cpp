#include "bandpass.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsbandpass {
namespace tool {
    
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
    params->progressCallback->cb = (rsReportProgressCallback_t) RSTool::showProgressCallback;
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

rsUIInterface* Bandpass::createUI()
{
    rsBandpassParameters *p = rsBandpassInitParameters();    
    rsBandpassBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rsbandpass::tool
