#include "rstimecourse.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rstimecourse {

RSTimecourse::RSTimecourse()
{}

void RSTimecourse::registerPlugin()
{
    RSTool::registerTool(createTimecourseToolRegistration());
    RSTool::registerXSDExtension(createTimecourseToolXSDExtension());
}

rsToolRegistration* RSTimecourse::createTimecourseToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->createTool = (rsToolToolCreator)RSTimecourse::createTimecourseTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSTimecourse::createTimecourseTask;
    return toolRegistration;
}

rsXSDExtension* RSTimecourse::createTimecourseToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rstimecourse.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSTimecourse::getName()
{
    return "Extract Timecourse";
}

const char* RSTimecourse::getCode()
{
    return "rstimecourse";
}

const char* RSTimecourse::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSTimecourse::createTimecourseTool()
{
    return (RSTool*) new tool::Timecourse();
}

RSTask* RSTimecourse::createTimecourseTask()
{
    return (RSTask*) new RSTask("rstimecourse", "Extract Timecourse");
}

}}}} // namespace rstools::batch::plugins::rstimecourse

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rstimecourse::RSTimecourse();
}
