#include "correlation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rscorrelation {
namespace tool {
    
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
    params->progressCallback->cb = (rsReportProgressCallback_t) RSTool::showProgressCallback;
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

rsUIInterface* Correlation::createUI()
{
    rsCorrelationParameters *p = rsCorrelationInitParameters();    
    rsCorrelationBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::correlation::tool
