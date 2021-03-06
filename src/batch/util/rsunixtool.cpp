#include <sys/stat.h>
#include "rsunixtool.hpp"
#include "rsunix.hpp"
#include <omp.h>
#include <nifti1_io.h>
#include <sys/types.h> 
#include <sys/wait.h>

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
    executionSuccessful = false;

    if (!this->_setupTempDir() || !this->_prepareStream() || !this->_prepareRun()) {
        return;
    }

    // Only fork when necessary (i.e. if we need to intercept the nifti output using streams)
    if (streamName == NULL) {
        // TODO: handle status properly
        const int returnStatus = system(executionCmd);
        executionSuccessful = returnStatus == 0;
    } else {

        // tokenize command
        const size_t argCount = 100;
        size_t i=1;
        char *args[argCount];
        char *cmd = rsString(executionCmd);
        if (!(args[0] = strtok(cmd, " \t"))) {
            fprintf(stderr, "Command could not be executed (invalid)\n");
            return;
        }
        while ((args[i] = strtok(NULL, " \t")) && i < argCount) {
            ++i;
        }
        args[i] = NULL;

        // create stream for intercepting the program's output
        FILE *stream = NULL;
        stream = fopen(streamName, "r+");

        // create thread that monitors the executing thread
        pid_t monitoringThreadPid;
        monitoringThreadPid = fork();

        if (monitoringThreadPid == 0) {

            // create execution thread
            pid_t executionThreadPid;
            executionThreadPid = fork();

            if (executionThreadPid == 0) {
                // child thread: execute the command
                execvp(args[0], args);
                fprintf(stderr, "Error while executing the command\n");
                exit(EXIT_FAILURE);
            }

            // monitoring thread: wait till the execution thread finishes
            int executionThreadStatus;
            waitpid(executionThreadPid, &executionThreadStatus, WUNTRACED);

            // ensure that the  main thread finishes in the case of an error in the execution
            // to do so, simply open and close the stream once more, so that the main thread
            // who's waiting for input on the stream receives a EOF
            if (stream != NULL) {
                FILE *writeStream = fopen(streamName, "w");
                if (writeStream) {
                    fclose(writeStream);
                }
            }

            // TODO: handle status properly
            const int executionThreadStatusCode = WEXITSTATUS(executionThreadStatus);
            if (WIFEXITED(executionThreadStatus) && executionThreadStatusCode == EXIT_SUCCESS) {
                exit(EXIT_SUCCESS);
            } else {
                exit(EXIT_FAILURE);
            }
        }

        // main thread: handle stream input
        if (stream != NULL) {
            executionSuccessful = this->_interceptOutput(stream);
        }

        int monitoringChildStatus;
        waitpid(monitoringThreadPid, &monitoringChildStatus, WUNTRACED);
        const int monitoringChildStatusCode = WEXITSTATUS(monitoringChildStatus);
        executionSuccessful &= WIFEXITED(monitoringChildStatus);
        executionSuccessful &= monitoringChildStatusCode == EXIT_SUCCESS;

        // close stream
        fclose(stream);
    }

    if (executionSuccessful) {
        this->_moveOutputIfNecessary();
    }
    this->_finalizeRun();
}

bool RSUnixTool::_interceptOutput(FILE *stream)
{
    // read in header
    const size_t headerSize = sizeof(nifti_1_header);
    nifti_1_header* header = (nifti_1_header*)rsMalloc(headerSize);

    if (fread(header, sizeof(nifti_1_header), 1L, stream) < 1) {
        fprintf(stderr, "Failed reading header information of the output nifti.\n");
        return false;
    }

    // read in nifti extender and extensions
    const size_t headerBytesLeft = header->vox_offset - headerSize;
    if (headerBytesLeft > 0) {
        char *buffer = (char*)rsMalloc(headerBytesLeft);
        const size_t remainingBytesRead = fread(buffer, sizeof(char), headerBytesLeft, stream);
        if (remainingBytesRead < headerBytesLeft) {
            fprintf(stderr, "Failed reading header extensions of the output nifti.\n");
            return false;
        }
        rsFree(buffer);
    }

    // replace extensions in the header we've received with the ones from the input
    rsAddCommentToNiftiHeader(inputNifti, executionCmdPrint);

    // calculate vox offset
    size_t voxOffset = headerSize + 4L; // start with 348 bytes for the header and 4 bytes for the nifti extender
    for (short i=0; i<inputNifti->num_ext; i++) {
        nifti1_extension e = inputNifti->ext_list[i];
        voxOffset += e.esize;
    }
    header->vox_offset = voxOffset;

    // write out header
    FILE *targetStream = fopen(streamTarget, "w+");

    if (targetStream < 0) {
        fprintf(stderr, "Failed opening the output nifti for writing (%s).\n", streamTarget);
        return false;
    }

    if (fwrite(header, headerSize, 1L, targetStream) < 1) {
        fclose(targetStream);
        fprintf(stderr, "Failed writing the output nifti (%s).\n", streamTarget);
        return false;
    }

    // write out nifti extender
    char niftiExtender[4] = {1, 0, 0, 0}; // will always be 1,0,0,0 as we'll always have at least one extension
    if (fwrite(niftiExtender, 4, 1L, targetStream) < 1) {
        fclose(targetStream);
        return false;
    }

    // write out extensions
    for (short i=0; i<inputNifti->num_ext; i++) {
        nifti1_extension e = inputNifti->ext_list[i];

        // write extension size
        if (fwrite(&e.esize, sizeof(int), 1L, targetStream) < 1) {
            fclose(targetStream);
            fprintf(stderr, "Failed writing the output nifti (%s). E1000\n", streamTarget);
            return false;
        }

        // write extension code
        if (fwrite(&e.ecode, sizeof(int), 1L, targetStream) < 1) {
            fclose(targetStream);
            fprintf(stderr, "Failed writing the output nifti (%s). E1001\n", streamTarget);
            return false;
        }

        // write extension size
        if (fwrite(e.edata, e.esize - 8L, 1L, targetStream) < 1) {
            fclose(targetStream);
            fprintf(stderr, "Failed writing the output nifti (%s). E1002\n", streamTarget);
            return false;
        }
    }

    // TODO: use dup2 instead of manually redirecting the stream
    // write data to the target destination by redirecting the input stream to the target stream
    //if (dup2(fileno(targetStream), fileno(stream)) < 0) {
    //    fclose(targetStream);
    //    return false;
    //}
    //
    //fflush(stream);

    // write out volume data
    const short *dims = &header->dim[1];
    const size_t expectedSize = rsGetBufferSize(dims[0], dims[1], dims[2], dims[3], header->datatype);
    const size_t blockSize = expectedSize / dims[3]; // read in one volume at a time
    void *buffer = rsMalloc(blockSize);
    size_t remainingBytes = expectedSize;
    rsFree(header);

    while (!feof(stream) && remainingBytes > 0) {
        const size_t actualBlockSize = min(remainingBytes, blockSize);
        const size_t nReadBytes = fread(buffer, 1L, actualBlockSize, stream);
        const size_t nWrittenBytes = fwrite(buffer, 1L, actualBlockSize, targetStream);
        if (nReadBytes < actualBlockSize || nWrittenBytes < actualBlockSize) {
            break;
        }
        remainingBytes -= nWrittenBytes;
    }

    rsFree(buffer);
    fclose(targetStream);

    if (remainingBytes < 1) {
        return true;
    } else {
        fprintf(stderr, "Failed writing the output nifti (%s). E1003\n", streamTarget);
        return false;
    }
}

bool RSUnixTool::_setupTempDir()
{
    // create temporary in which everything gets executed
    char* tmpDirNameTpl = (char*)malloc(sizeof(char)*255);
    sprintf(tmpDirNameTpl, "%s", "/tmp/rsbatch.unix.cmd-XXXXXXX");
    tmpDirPath = mkdtemp(tmpDirNameTpl);

    if (tmpDirPath == NULL) {
        fprintf(stderr, "Could not create a temporary directory for the execution of a command line script.\n");
        return false;
    }

    getUnixTask()->setTempDirectoryPath(rsString(tmpDirPath));

    return true;
}

bool RSUnixTool::_prepareRun()
{
    // place a wrapper script that calls the command to be executed and redirects its output
    char* wrapperScriptName = rsStringConcat(tmpDirPath, "/cmdWrapper.sh", NULL);
    ofstream wrapperScript;
    wrapperScript.open(wrapperScriptName);
    wrapperScript << "#!/bin/bash" << std::endl;
    wrapperScript << tmpDirPath << "/cmd.sh 2> " << tmpDirPath <<  "/error.log 1> " << tmpDirPath <<  "/output.log" << std::endl;
    wrapperScript << "returnCode=\"${PIPESTATUS[0]}\"" << std::endl;
    wrapperScript << "cat " << tmpDirPath <<  "/error.log >&2" << std::endl;
    wrapperScript << "exit $returnCode" << std::endl;
    wrapperScript.close();

    // place a bash script containing the command to be executed in the directory
    char* scriptName = rsStringConcat(tmpDirPath, "/cmd.sh", NULL);
    ofstream script;
    script.open(scriptName);
    script << "#!/bin/bash" << std::endl;
    script << this->getUnixTask()->getCmd(true) << std::endl;
    script << "exit 0" << std::endl;
    script.close();

    // fix permissions
    chmod(wrapperScriptName, S_IXUSR | S_IWUSR | S_IRUSR | S_IRGRP);
    chmod(scriptName, S_IXUSR | S_IWUSR | S_IRUSR | S_IRGRP);

    // prepare statement to call it with
    executionCmd = wrapperScriptName;
    rsFree(scriptName);

    executionCmdPrint = this->getUnixTask()->getCmd(false);

    return true;
}

bool RSUnixTool::_prepareStream()
{
    // create stream to intercept the nifti that is created and alter its header before it is written out
    // if requested by the individual tool extending this class
    streamName = NULL;
    return true;
}

bool RSUnixTool::_createStream(const char* path)
{
    streamName = rsStringConcat(tmpDirPath, "/stream.nii", NULL);
    if (mkfifo(path, 0666) != 0) {
        fprintf(stderr, "Could not create a temporary directory for the execution of a command line script.\n");
        return false;
    }

    return true;
}

void RSUnixTool::_moveOutputIfNecessary()
{
    // if required a tool extending this class can move some additionally created files
    // from the temporary created folder by overwriting this method
}

void RSUnixTool::_finalizeRun()
{
    // if it was executed correctly destroy the temporary working directory, etc.
    if ( executionSuccessful ) {
        char *rmCommand = rsStringConcat("rm -rf ", tmpDirPath, NULL);
        system(rmCommand);
        rsFree(rmCommand);
    } else { // otherwise keep the dir so that the user can debug it
        fprintf(stderr, "Error while executing shell task. For more information see '%s'.\n", tmpDirPath);
    }

    rsFree(executionCmd);
    rsFree(executionCmdPrint);
    rsFree(streamName);
    rsFree(tmpDirPath);
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
    
    fprintf(stream, "Cmd:\n%s\n", getUnixTask()->getCmd(false));
    fprintf(stream, "\n");
}

}}} // namespace rstools::batch::util
