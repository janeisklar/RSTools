#include "rsroi.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsroi {

RSRoi::RSRoi()
{}

void RSRoi::registerPlugin()
{
    RSTool::registerTool(createRoiToolRegistration());
    RSTool::registerXSDExtension(createRoiToolXSDExtension());
}

rsToolRegistration* RSRoi::createRoiToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->createTool = (rsToolToolCreator)RSRoi::createRoiTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSRoi::createRoiTask;
    return toolRegistration;
}

rsXSDExtension* RSRoi::createRoiToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rsroi.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSRoi::getName()
{
    return "Create ROI";
}

const char* RSRoi::getCode()
{
    return "rsroi";
}

const char* RSRoi::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSRoi::createRoiTool()
{
    return (RSTool*) new tool::Roi();
}

RSTask* RSRoi::createRoiTask()
{
    return (RSTask*) new RSTask("rsroi", "Create ROI");
}

}}}} // namespace rstools::batch::plugins::rsroi

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rsroi::RSRoi();
}
