#include "applytransformation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsapplytransformation {
namespace tool {
    
void ApplyTransformation::_parseParams(int argc, char * argv[])
{
    params = rsApplyTransformationParseParams(argc, argv);
}
    
void ApplyTransformation::_init()
{
    params->threads = this->threads;
    rsApplyTransformationInit(params);
}

void ApplyTransformation::_run()
{
    params->progressCallback = (rsReportProgressCallback*)rsMalloc(sizeof(rsReportProgressCallback));
    params->progressCallback->cb = (rsReportProgressCallback_t) RSTool::showProgressCallback;
    params->progressCallback->data = (void*)oc;

    rsApplyTransformationRun(params);
    
    rsFree(params->progressCallback);
}

void ApplyTransformation::destroy()
{
    rsApplyTransformationDestroy(params);
}

bool ApplyTransformation::isEverythingFine()
{
    return params != NULL && params->parametersValid;
}

rsUIInterface* ApplyTransformation::createUI()
{
    rsApplyTransformationParameters *p = rsApplyTransformationInitParameters();
    rsApplyTransformationBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rsapplytransformation::tool
