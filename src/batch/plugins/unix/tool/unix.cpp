#include "unix.hpp"
#include "utils/rsstring.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>

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
    this->executionSuccessful = false;
    
    // create a temporay directory where we store the script and log messages
    char* tmpDirNameTpl = (char*)malloc(sizeof(char)*255);
    sprintf(tmpDirNameTpl, "%s", "/tmp/rsbatch.unix.cmd-XXXXXXX");
    char* tmpDirName = mkdtemp(tmpDirNameTpl);

    if ( tmpDirName == NULL ) {
        throw runtime_error("Could not create a temporary directory for the execution of a command line script.");
    }
    
    // place a bash script containing the command to be executed in the directory
    char* scriptName = rsStringConcat(tmpDirName, "/cmd.sh", NULL);
    ofstream script;
    script.open(scriptName);
    script << this->getUnixTask()->getCmd();
    script.close();
    
    // prepare statement to call it with
    char* executionCmd = rsStringConcat("/bin/bash ", scriptName, " > ", tmpDirName, "/output.log", NULL);
    
    // if it was executed correctly destroy the temporary working directory
    if ( system(executionCmd) == 0 ) {
        this->executionSuccessful = true;
        
        char *rmCommand = rsStringConcat("rm -rf ", tmpDirName, NULL);
        system(rmCommand);
        rsFree(rmCommand);
    } else { // otherwise keep the dir so that the user can debug it
        fprintf(stderr, "Error while executing shell task. For more information see '%s'.\n", tmpDirName);
    }
    
    rsFree(executionCmd);
    rsFree(scriptName);
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
