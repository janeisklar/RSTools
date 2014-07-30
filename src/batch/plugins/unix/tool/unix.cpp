#include "unix.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {
namespace tool {
    
void Unix::_parseParams(int argc, char * argv[])
{   
    this->executionSuccessful = true;
}
    
void Unix::_init()
{}

void Unix::_run()
{
    if ( system(this->getUnixTask()->getCmd()) != 0 ) {
        this->executionSuccessful = false;
    }
}

void Unix::destroy()
{}

bool Unix::isEverythingFine()
{
    return this->executionSuccessful;
}

rstools::batch::plugins::unix::task::Unix* Unix::getUnixTask()
{
    return (rstools::batch::plugins::unix::task::Unix*)this->getTask();
}

rsUIInterface* Unix::createUI()
{
    rsUIOption *o;
    rsUIInterface* interface = rsUINewInterface();
    interface->description   = "Execute Unix Command";
    
    o = rsUINewOption();
    o->name                = "command";
    o->type                = G_OPTION_ARG_STRING;
    o->cli_description     = "the unix command that is to be executed";
    o->cli_arg_description = "<unix cmd>";
    rsUIAddOption(interface, o);
    
    return interface;
}

}}}}} // namespace rstools::batch::plugins::unix::tool
