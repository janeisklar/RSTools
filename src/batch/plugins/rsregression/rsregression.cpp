#include "rsregression.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsregression {

RSRegression::RSRegression()
{
    
}
    
void RSRegression::registerPlugin()
{
    RSTool::registerTool(createRegressionToolRegistration());
    RSTool::registerXSDExtension(createRegressionToolXSDExtension());
}

rsToolRegistration* RSRegression::createRegressionToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->createTool = (rsToolToolCreator)RSRegression::createRegressionTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSRegression::createRegressionTask;
    return toolRegistration;
}
rsXSDExtension* RSRegression::createRegressionToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rsregression.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSRegression::getName()
{
    return "Regression";
}

const char* RSRegression::getCode()
{
    return "rsregression";
}

const char* RSRegression::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSRegression::createRegressionTool()
{
    return (RSTool*) new tool::Regression();
}

RSTask* RSRegression::createRegressionTask()
{
    return (RSTask*) new RSTask("rsregression", "Regression");
}

}}}} // namespace rstools::batch::plugins::rsregression

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rsregression::RSRegression();
}
