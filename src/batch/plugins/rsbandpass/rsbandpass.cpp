#include "rsbandpass.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsbandpass {

RSBandpass::RSBandpass()
{}

void RSBandpass::registerPlugin()
{
    RSTool::registerTool(createBandpassToolRegistration());
    RSTool::registerXSDExtension(createBandpassToolXSDExtension());
}

rsToolRegistration* RSBandpass::createBandpassToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->createTool = (rsToolToolCreator)RSBandpass::createBandpassTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSBandpass::createBandpassTask;
    return toolRegistration;
}
rsXSDExtension* RSBandpass::createBandpassToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rsbandpass.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSBandpass::getName()
{
    return "Bandpass Filter";
}

const char* RSBandpass::getCode()
{
    return "rsbandpass";
}

const char* RSBandpass::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSBandpass::createBandpassTool()
{
    return (RSTool*) new tool::Bandpass();
}

RSTask* RSBandpass::createBandpassTask()
{
    return (RSTask*) new RSTask("rsbandpass", "Bandpass Filter");
}

}}}} // namespace rstools::batch::plugins::rsbandpass

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rsbandpass::RSBandpass();
}
