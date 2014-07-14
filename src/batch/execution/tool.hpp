#ifndef rstools_rsbatch_execution_tool_h
#define rstools_rsbatch_execution_tool_h

#include <iostream>
#include "src/batch/outputCatcher.hpp"
#include "src/batch/rstask.hpp"
#include <sys/ioctl.h>
#include <math.h>
#include <stdexcept>
#include "src/rscommon.h"

namespace rstools {
namespace batch {
namespace execution {

class Tool {

    public:
        static Tool* factory(const short taskCode);
        
        void parseParams(int argc, char** argv);
        void init();
        void run();
        virtual void destroy() = 0;
        virtual bool isEverythingFine() = 0;
        
        char const * getOutput();
        char **getCallString(int *argc);
        RSTask *getTask();
        void setTask(RSTask *task);
        void setThreads(int threads);
        
        static void showProgressCallback(rsReportProgressEvent *event, void *userdata);
        static void printProgressBar(FILE* stream, double percentage, int run, char* description);
        
        static const short tools[];
        
    protected:
        virtual void _parseParams(int, char**) = 0;
        virtual void _init() = 0;
        virtual void _run() = 0;
        std::string _message;
        int argc;
        char** argv;
        RSTask *task;
        int threads;
        rstools::batch::OutputCatcher *oc;
};

}}} // namespace rstools::batch::execution

#endif
