#include "smoothing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rssmoothing {
namespace tool {
    
void Smoothing::_parseParams(int argc, char * argv[])
{
    params = rsSmoothingParseParams(argc, argv);
}
    
void Smoothing::_init()
{
    params->threads = this->threads;    
    rsSmoothingInit(params);
}

void Smoothing::_run()
{
    params->progressCallback = (rsReportProgressCallback*)rsMalloc(sizeof(rsReportProgressCallback));
    params->progressCallback->cb = (rsReportProgressCallback_t) RSTool::showProgressCallback;
    params->progressCallback->data = (void*)oc;
        
    rsSmoothingRun(params);
    
    rsFree(params->progressCallback);
}

void Smoothing::destroy()
{
    rsSmoothingDestroy(params);
}

bool Smoothing::isEverythingFine()
{
    return params != NULL && params->parametersValid;
}

rsUIInterface* Smoothing::createUI()
{
    rsSmoothingParameters *p = rsSmoothingInitParameters();    
    rsSmoothingBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rssmoothing::tool
