#include "regression.hpp"

namespace rstools {
namespace batch {
namespace execution {
    
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
    params->progressCallback->cb = (rsReportProgressCallback_t) Tool::showProgressCallback;
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

}}} // namespace rstools::batch::execution
