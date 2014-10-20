#include "rsunixtool.hpp"
#include "rsunix.hpp"

namespace rstools {
namespace batch {
namespace util {


void RSUnixTool::_parseParams(int argc, char * argv[])
{   
    this->executionSuccessful = true;
}

bool RSUnixTool::isEverythingFine()
{
    return this->executionSuccessful;
}

RSUnixTask* RSUnixTool::getUnixTask()
{
    return (RSUnixTask*)this->getTask();
}

void RSUnixTool::_run()
{
    this->executionSuccessful = rsExecuteUnixCommand(this->getUnixTask()->getCmd());
}

void RSUnixTool::printCallString(FILE *stream)
{
    int argc;
    char **argv = getCallString(&argc);

    fprintf(stream, "Tool:\n %s\n\nParams:\n", getTask()->getName());
    for ( int i=1; i<argc; i++ ) {
        fprintf(stream, "  %s\n", argv[i]);
    }
    
    fprintf(stream, "\n");
    
    fprintf(stream, "Cmd:\n%s\n", getUnixTask()->getCmd());
    fprintf(stream, "\n");
}

}}} // namespace rstools::batch::util
