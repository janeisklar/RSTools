#ifndef rstools_rsbatch_util_tool_h
#define rstools_rsbatch_util_tool_h

#include <iostream>
#include <sys/ioctl.h>
#include <math.h>
#include <stdexcept>
#include <vector>
#include "rscommon.h"
#include "outputCatcher.hpp"
#include "rstask.hpp"
#include "utils/rsui.h"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

class RSTask;
class RSTool;

typedef RSTool* (*rsToolToolCreator) (void);
typedef RSTask* (*rsToolTaskCreator) (void);

typedef struct {
    
    // name of the tool that will be used
    // in the GUI/batch processor
    const char* name;
    
    // internal code that refers to the
    // tool
    const char* code;
    
    // callback function that creates a
    // new instance of the tool
    rsToolToolCreator createTool;
    
    // callback function that creates a
    // new instance of the task for the
    // tool
    rsToolTaskCreator createTask;
    
} rsToolRegistration;

typedef struct {
    
    // name of the element in the job-file
    // should be the same as the code for 
    // the task/tool
    const char* name;
    
    // file that will be merged into the
    // main job-XSD, which should contain
    // an XSD type definition, which should
    // be derived from the RSTask type as
    // defined in the main job-XSD
    const char* file;
    
    // name of the XSD-type that will be
    // used in the merged XSD for referencing
    // the new XML node type
    const char* type;
    
} rsXSDExtension;

class RSTool {

    public:
        static rsToolRegistration* findRegistration(const char* code);
        static RSTool* toolFactory(const char* code);
        static void registerTool(rsToolRegistration* registration);
        static void registerXSDExtension(rsXSDExtension* extension);
        
        virtual void parseParams(int argc, char** argv);
        virtual void init();
        virtual void run();
        virtual void destroy() = 0;
        virtual bool isEverythingFine() = 0;
        virtual rsUIInterface* createUI() = 0;
        
        virtual char const * getOutput();
        virtual char **getCallString(int *argc);
        virtual RSTask *getTask();
        virtual void setTask(RSTask *task);
        virtual void setThreads(int threads);
        
        virtual void printCallString(FILE *stream);
        
        static void showProgressCallback(rsReportProgressEvent *event, void *userdata);
        static void printProgressBar(FILE* stream, double percentage, int run, char* description);
        
        static vector<const char*> getTools();
        static vector<rsXSDExtension*> getXSDExtensions();
        
    protected:
        virtual void _parseParams(int, char**) = 0;
        virtual void _init() = 0;
        virtual void _run() = 0;
        
        string _message;
        int argc;
        char** argv;
        RSTask *task;
        int threads;
        OutputCatcher *oc;
        static vector<rsToolRegistration*> tools;
        static vector<rsXSDExtension*> xsdExtensions;
};

}}} // namespace rstools::batch::util

#endif
