#include "rssmoothing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rssmoothing {

RSSmoothing::RSSmoothing()
{}

void RSSmoothing::registerPlugin()
{
    RSTool::registerTool(createSmoothingToolRegistration());
    RSTool::registerXSDExtension(createSmoothingToolXSDExtension());
}

rsToolRegistration* RSSmoothing::createSmoothingToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->category   = "Spatial";
    toolRegistration->createTool = (rsToolToolCreator)RSSmoothing::createSmoothingTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSSmoothing::createSmoothingTask;
    return toolRegistration;
}
rsXSDExtension* RSSmoothing::createSmoothingToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rssmoothing.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSSmoothing::getName()
{
    return "Smoothing";
}

const char* RSSmoothing::getCode()
{
    return "rssmoothing";
}

const char* RSSmoothing::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSSmoothing::createSmoothingTool()
{
    return (RSTool*) new tool::Smoothing();
}

RSTask* RSSmoothing::createSmoothingTask()
{
    return (RSTask*) new RSTask("rssmoothing", "Smoothing");
}

}}}} // namespace rstools::batch::plugins::rssmoothing

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rssmoothing::RSSmoothing();
}
