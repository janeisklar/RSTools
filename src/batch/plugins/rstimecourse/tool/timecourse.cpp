#include "timecourse.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rstimecourse {
namespace tool {
    
void Timecourse::_parseParams(int argc, char * argv[])
{
    params = rsTimecourseParseParams(argc, argv);
}
    
void Timecourse::_init()
{
    params->threads = this->threads;
    rsTimecourseInit(params);
}

void Timecourse::_run()
{
    rsTimecourseRun(params);
}

void Timecourse::destroy()
{
    rsTimecourseDestroy(params);
}

bool Timecourse::isEverythingFine()
{
    return params != NULL && params->parametersValid;
}

rsUIInterface* Timecourse::createUI()
{
    rsTimecourseParameters *p = rsTimecourseInitParameters();    
    rsTimecourseBuildInterface(p);
    
    return p->interface;
}

}}}}} // namespace rstools::batch::plugins::rstimecourse::tool
