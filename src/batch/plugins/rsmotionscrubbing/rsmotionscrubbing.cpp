#include "rsmotionscrubbing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsmotionscrubbing {

RSMotionScrubbing::RSMotionScrubbing()
{}

void RSMotionScrubbing::registerPlugin()
{
    RSTool::registerTool(createMotionScrubbingToolRegistration());
    RSTool::registerXSDExtension(createMotionScrubbingToolXSDExtension());
}

rsToolRegistration* RSMotionScrubbing::createMotionScrubbingToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->category   = "Artifacts / Motion";
    toolRegistration->createTool = (rsToolToolCreator)RSMotionScrubbing::createMotionScrubbingTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSMotionScrubbing::createMotionScrubbingTask;
    return toolRegistration;
}

rsXSDExtension* RSMotionScrubbing::createMotionScrubbingToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rsmotionscrubbing.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSMotionScrubbing::getName()
{
    return "Motion Scrubbing";
}

const char* RSMotionScrubbing::getCode()
{
    return "rsmotionscrubbing";
}

const char* RSMotionScrubbing::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSMotionScrubbing::createMotionScrubbingTool()
{
    return (RSTool*) new tool::MotionScrubbing();
}

RSTask* RSMotionScrubbing::createMotionScrubbingTask()
{
    return (RSTask*) new RSTask("rsmotionscrubbing", "Motion Scrubbing");
}

}}}} // namespace rstools::batch::plugins::rsmotionscrubbing

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rsmotionscrubbing::RSMotionScrubbing();
}
