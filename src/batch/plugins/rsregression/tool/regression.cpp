#include "regression.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsregression {
namespace tool {
    
void Regression::_parseParams(int argc, char * argv[])
{
    params = rsRegressionParseParams(argc, argv);
}
    
void Regression::_init()
{
    params->threads = this->threads;
    rsRegressionInit(params);
}

void Regression::_run()
{
    params->progressCallback = (rsReportProgressCallback*)rsMalloc(sizeof(rsReportProgressCallback));
    params->progressCallback->cb = (rsReportProgressCallback_t) RSTool::showProgressCallback;
    params->progressCallback->data = (void*)oc;
        
    rsRegressionRun(params);
    
    rsFree(params->progressCallback);
}

void Regression::destroy()
{
    rsRegressionDestroy(params);
}

bool Regression::isEverythingFine()
{
    return params != NULL && params->parametersValid;
}

rsUIInterface* Regression::createUI()
{
    rsRegressionParameters *p = rsRegressionInitParameters();    
    rsRegressionBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rsregression::tool
