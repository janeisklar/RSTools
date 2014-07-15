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
    fprintf(stdout, "Hello world regression!! \n");
}

const char* RSRegression::getName()
{
    return "Regression";
}

const char* RSRegression::getVersion()
{
    return RSTOOLS_VERSION_LABEL;
}

}}}} // namespace rstools::batch::plugins::rsregression

Plugin* rsGetPlugin(void)
{
    return (Plugin*) new rsregression::RSRegression();
}
