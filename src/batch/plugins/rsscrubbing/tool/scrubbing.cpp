#include "scrubbing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsscrubbing {
namespace tool {
    
void Scrubbing::_parseParams(int argc, char * argv[])
{
    params = rsScrubbingParseParams(argc, argv);
}
    
void Scrubbing::_init()
{
    rsScrubbingInit(params);
}

void Scrubbing::_run()
{
    rsScrubbingRun(params);
}

void Scrubbing::destroy()
{
    rsScrubbingDestroy(params);
}

bool Scrubbing::isEverythingFine()
{
    return params != NULL && params->parametersValid;
}

rsUIInterface* Scrubbing::createUI()
{
    rsScrubbingParameters *p = rsScrubbingInitParameters();    
    rsScrubbingBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rsscrubbing::tool
