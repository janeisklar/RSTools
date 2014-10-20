#include "unix.hpp"
#include "utils/rsstring.h"
#include "batch/util/rsunix.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {
namespace tool {
    
void Unix::_init()
{}

void Unix::destroy()
{}

rsUIInterface* Unix::createUI()
{
    rsUIOption *o;
    rsUIInterface* interface = rsUINewInterface();
    interface->description   = rsString("Execute Unix Command");
    
    o = rsUINewOption();
    o->name                = rsString("command");
    o->type                = G_OPTION_ARG_STRING;
    o->cli_description     = rsString("the unix command that is to be executed");
    o->cli_arg_description = rsString("<unix cmd>");
    o->nLines              = 20;
    rsUIAddOption(interface, o);
    
    return interface;
}

void Unix::printCallString(FILE *stream)
{
    int argc;
    char **argv = getCallString(&argc);

    fprintf(stream, "Tool:\n %s\n\n", getTask()->getName());
    fprintf(stream, "Cmd:\n%s\n", getUnixTask()->getCmd());
    fprintf(stream, "\n");
}

}}}}} // namespace rstools::batch::plugins::unix::tool
