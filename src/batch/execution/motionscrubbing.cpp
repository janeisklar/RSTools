#include "motionscrubbing.hpp"

namespace rstools {
namespace batch {
namespace execution {
    
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

}}} // namespace rstools::batch::execution
