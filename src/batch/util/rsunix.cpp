#include "rsunix.hpp"
#include "utils/rsstring.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

using namespace std;

bool rsExecuteUnixCommand(const char *cmd, char *tmpDir)
{
    bool executionSuccessful = false;
    char *tmpDirName;

    if (tmpDir == NULL) {
        // create a temporay directory where we store the script and log messages
        char *tmpDirNameTpl = (char *) malloc(sizeof(char) * 255);
        sprintf(tmpDirNameTpl, "%s", "/tmp/rsbatch.unix.cmd-XXXXXXX");
        tmpDirName = mkdtemp(tmpDirNameTpl);

        if (tmpDirName == NULL) {
            throw runtime_error("Could not create a temporary directory for the execution of a command line script.");
        }
    } else {
        tmpDirName = tmpDir;
    }
    
    // place a bash script containing the command to be executed in the directory
    char* scriptName = rsStringConcat(tmpDirName, "/cmd.sh", NULL);
    ofstream script;
    script.open(scriptName);
    script << cmd;
    script.close();
    
    // prepare statement to call it with
    char* executionCmd = rsStringConcat("/bin/bash ", scriptName, " > ", tmpDirName, "/output.log", NULL);
    
    // if it was executed correctly destroy the temporary working directory
    if ( system(executionCmd) == 0 ) {
        executionSuccessful = true;
        
        char *rmCommand = rsStringConcat("rm -rf ", tmpDirName, NULL);
        system(rmCommand);
        rsFree(rmCommand);
    } else { // otherwise keep the dir so that the user can debug it
        fprintf(stderr, "Error while executing shell task. For more information see '%s'.\n", tmpDirName);
    }
    
    rsFree(executionCmd);
    rsFree(scriptName);
    
    return executionSuccessful;
}
