#include "orientation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsorientation {
namespace tool {
    
void Orientation::_parseParams(int argc, char * argv[])
{
    params = rsOrientationParseParams(argc, argv);
}
    
void Orientation::_init()
{
    rsOrientationInit(params);
}

void Orientation::_run()
{
    rsOrientationRun(params);
}

void Orientation::destroy()
{
    rsOrientationDestroy(params);
}

bool Orientation::isEverythingFine()
{
    return params != NULL && params->parametersValid;
}

rsUIInterface* Orientation::createUI()
{
    rsOrientationParameters *p = rsOrientationInitParameters();    
    rsOrientationBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rsorientation::tool
