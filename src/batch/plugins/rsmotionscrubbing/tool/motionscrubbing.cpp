#include "motionscrubbing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsmotionscrubbing {
namespace tool {
    
void MotionScrubbing::_parseParams(int argc, char * argv[])
{
    params = rsMotionScrubbingParseParams(argc, argv);
}
    
void MotionScrubbing::_init()
{
    params->threads = this->threads;    
    rsMotionScrubbingInit(params);
}

void MotionScrubbing::_run()
{
    rsMotionScrubbingRun(params);
}

void MotionScrubbing::destroy()
{
    rsMotionScrubbingDestroy(params);
}

bool MotionScrubbing::isEverythingFine()
{
    return params != NULL && params->parametersValid;
}

rsUIInterface* MotionScrubbing::createUI()
{
    rsMotionScrubbingParameters *p = rsMotionScrubbingInitParameters();    
    rsMotionScrubbingBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rsmotionscrubbing::tool
