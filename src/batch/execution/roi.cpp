#include "roi.hpp"

namespace rstools {
namespace batch {
namespace execution {
    
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

}}} // namespace rstools::batch::execution
