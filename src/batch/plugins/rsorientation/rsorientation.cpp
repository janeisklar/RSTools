#include "rsorientation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsorientation {

RSOrientation::RSOrientation()
{}

void RSOrientation::registerPlugin()
{
    RSTool::registerTool(createOrientationToolRegistration());
    RSTool::registerXSDExtension(createOrientationToolXSDExtension());
}

rsToolRegistration* RSOrientation::createOrientationToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->category   = "Spatial";
    toolRegistration->createTool = (rsToolToolCreator)RSOrientation::createOrientationTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSOrientation::createOrientationTask;
    return toolRegistration;
}
rsXSDExtension* RSOrientation::createOrientationToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rsorientation.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSOrientation::getName()
{
    return "Orientation";
}

const char* RSOrientation::getCode()
{
    return "rsorientation";
}

const char* RSOrientation::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSOrientation::createOrientationTool()
{
    return (RSTool*) new tool::Orientation();
}

RSTask* RSOrientation::createOrientationTask()
{
    return (RSTask*) new RSTask("rsorientation", "Orientation");
}

}}}} // namespace rstools::batch::plugins::rsorientation

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rsorientation::RSOrientation();
}
