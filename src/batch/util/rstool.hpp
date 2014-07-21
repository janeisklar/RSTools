#ifndef rstools_rsbatch_util_tool_h
#define rstools_rsbatch_util_tool_h

#include <iostream>
#include <sys/ioctl.h>
#include <math.h>
#include <stdexcept>
#include <vector>
#include "src/rscommon.h"
#include "outputCatcher.hpp"
#include "rstask.hpp"
#include "src/utils/rsui.h"

namespace rstools {
namespace batch {
namespace util {

class RSTool;

typedef RSTool* (*rsToolToolCreator) (void);
typedef RSTask* (*rsToolTaskCreator) (void);

typedef struct {
    const char* name;
    const char* code;
    rsToolToolCreator createTool;
    rsToolTaskCreator createTask;
} rsToolRegistration;

class RSTool {

    public:
        static rsToolRegistration* findRegistration(const char* code);
        static RSTool* toolFactory(const char* code);
        static void registerTool(rsToolRegistration* registration);
        
        void parseParams(int argc, char** argv);
        void init();
        void run();
        virtual void destroy() = 0;
        virtual bool isEverythingFine() = 0;
        virtual rsUIInterface* createUI() = 0;
        
        char const * getOutput();
        char **getCallString(int *argc);
        RSTask *getTask();
        void setTask(RSTask *task);
        void setThreads(int threads);
        
        void printCallString(FILE *stream);
        
        static void showProgressCallback(rsReportProgressEvent *event, void *userdata);
        static void printProgressBar(FILE* stream, double percentage, int run, char* description);
        
        static vector<const char*> getTools();
        
    protected:
        virtual void _parseParams(int, char**) = 0;
        virtual void _init() = 0;
        virtual void _run() = 0;
        std::string _message;
        int argc;
        char** argv;
        RSTask *task;
        int threads;
        OutputCatcher *oc;
        static vector<rsToolRegistration*> tools;
};

}}} // namespace rstools::batch::util

#endif
