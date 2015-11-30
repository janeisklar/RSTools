#include "rsapplytransformation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsapplytransformation {

RSApplyTransformation::RSApplyTransformation()
{}

void RSApplyTransformation::registerPlugin()
{
    RSTool::registerTool(createApplyTransformationToolRegistration());
    RSTool::registerXSDExtension(createApplyTransformationToolXSDExtension());
}

rsToolRegistration* RSApplyTransformation::createApplyTransformationToolRegistration()
{
    rsToolRegistration* toolRegistration = (rsToolRegistration*)malloc(sizeof(rsToolRegistration));
    toolRegistration->name       = getName();
    toolRegistration->code       = getCode();
    toolRegistration->category   = "Spatial";
    toolRegistration->createTool = (rsToolToolCreator)RSApplyTransformation::createApplyTransformationTool;
    toolRegistration->createTask = (rsToolTaskCreator)RSApplyTransformation::createApplyTransformationTask;
    return toolRegistration;
}
rsXSDExtension* RSApplyTransformation::createApplyTransformationToolXSDExtension()
{
    rsXSDExtension* toolExtension = (rsXSDExtension*)malloc(sizeof(rsXSDExtension));
    toolExtension->name           = getCode();
    toolExtension->file           = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/plugins/rsapplytransformation.xsdext";
    toolExtension->type           = getCode();
    return toolExtension;
}

const char* RSApplyTransformation::getName()
{
    return "Unified Transformation";
}

const char* RSApplyTransformation::getCode()
{
    return "rsapplytransformation";
}

const char* RSApplyTransformation::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

RSTool* RSApplyTransformation::createApplyTransformationTool()
{
    return (RSTool*) new tool::ApplyTransformation();
}

RSTask* RSApplyTransformation::createApplyTransformationTask()
{
    return (RSTask*) new RSTask("rsapplytransformation", "Unified Transformation");
}

}}}} // namespace rstools::batch::plugins::rsapplytransformation

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rstools::batch::plugins::rsapplytransformation::RSApplyTransformation();
}
