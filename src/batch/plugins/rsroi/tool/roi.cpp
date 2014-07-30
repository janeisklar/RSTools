#include "roi.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsroi {
namespace tool {
    
void Roi::_parseParams(int argc, char * argv[])
{
    params = rsRoiParseParams(argc, argv);
}
    
void Roi::_init()
{
    rsRoiInit(params);
}

void Roi::_run()
{
    rsRoiRun(params);
}

void Roi::destroy()
{
    rsRoiDestroy(params);
}

bool Roi::isEverythingFine()
{
    return params != NULL && params->parametersValid;
}

rsUIInterface* Roi::createUI()
{
    rsRoiParameters *p = rsRoiInitParameters();    
    rsRoiBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rsroi::tool
