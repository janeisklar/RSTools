#include "rscorrelation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rscorrelation {

RSCorrelation::RSCorrelation()
{
    
}
    
void RSCorrelation::registerPlugin()
{
    RSTool::registerTool(createCorrelationToolRegistration());
    RSTool::registerXSDExtension(createCorrelationToolXSDExtension());
}

rsToolRegistration* RSCorrelation::createCorrelationToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->category   = "Temporal";
    toolRegistration->createTool = (rsToolToolCreator)RSCorrelation::createCorrelationTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSCorrelation::createCorrelationTask;
    return toolRegistration;
}

rsXSDExtension* RSCorrelation::createCorrelationToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rscorrelation.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSCorrelation::getName()
{
    return "Correlation";
}

const char* RSCorrelation::getCode()
{
    return "rscorrelation";
}

const char* RSCorrelation::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSCorrelation::createCorrelationTool()
{
    return (RSTool*) new tool::Correlation();
}

RSTask* RSCorrelation::createCorrelationTask()
{
    return (RSTask*) new RSTask("rscorrelation", "Correlation");
}

}}}} // namespace rstools::batch::plugins::rscorrelation

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rscorrelation::RSCorrelation();
}
