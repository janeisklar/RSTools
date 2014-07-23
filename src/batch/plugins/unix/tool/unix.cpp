#include "unix.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {
namespace tool {
    
void Unix::_parseParams(int argc, char * argv[])
{
    fprintf(stdout, "Cmd: %s\n", argv[1]);
    
    #pragma message "TODO: Implement Unix::_parseParams() in " __FILE__
    
}
    
void Unix::_init()
{}

void Unix::_run()
{
    #pragma message "TODO: Implement Unix::_run() in " __FILE__
}

void Unix::destroy()
{}

bool Unix::isEverythingFine()
{
    return true;
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
