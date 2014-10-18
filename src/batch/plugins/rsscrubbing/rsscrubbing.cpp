#include "rsscrubbing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsscrubbing {

RSScrubbing::RSScrubbing()
{}

void RSScrubbing::registerPlugin()
{
    RSTool::registerTool(createScrubbingToolRegistration());
    RSTool::registerXSDExtension(createScrubbingToolXSDExtension());
}

rsToolRegistration* RSScrubbing::createScrubbingToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->category   = "Artifacts / Motion";
    toolRegistration->createTool = (rsToolToolCreator)RSScrubbing::createScrubbingTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSScrubbing::createScrubbingTask;
    return toolRegistration;
}

rsXSDExtension* RSScrubbing::createScrubbingToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rsscrubbing.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSScrubbing::getName()
{
    return "Apply Scrubbing";
}

const char* RSScrubbing::getCode()
{
    return "rsscrubbing";
}

const char* RSScrubbing::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSScrubbing::createScrubbingTool()
{
    return (RSTool*) new tool::Scrubbing();
}

RSTask* RSScrubbing::createScrubbingTask()
{
    return (RSTask*) new RSTask("rsscrubbing", "Apply Scrubbing");
}

}}}} // namespace rstools::batch::plugins::rsscrubbing

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rsscrubbing::RSScrubbing();
}
