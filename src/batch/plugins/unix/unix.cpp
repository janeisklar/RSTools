#include "unix.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {

Unix::Unix()
{}

void Unix::registerPlugin()
{
    RSTool::registerTool(createUnixToolRegistration());
    RSTool::registerXSDExtension(createUnixToolXSDExtension());
}

rsToolRegistration* Unix::createUnixToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->category   = "Util";
    toolRegistration->createTool = (rsToolToolCreator)Unix::createUnixTool;
    toolRegistration->createTask = (rsToolTaskCreator)Unix::createUnixTask;
    return toolRegistration;
}
rsXSDExtension* Unix::createUnixToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/unix.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* Unix::getName()
{
    return "Execute Unix Command";
}

const char* Unix::getCode()
{
    return "unix";
}

const char* Unix::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* Unix::createUnixTool()
{
    return (RSTool*) new tool::Unix();
}

RSTask* Unix::createUnixTask()
{
    return (RSTask*) new task::Unix("unix", "Execute Unix Command");
}

}}}} // namespace rstools::batch::plugins::unix

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::unix::Unix();
}
