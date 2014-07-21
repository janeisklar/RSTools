#ifndef rstools_rsbatch_util_rstask_hpp
#define rstools_rsbatch_util_rstask_hpp

#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <locale>
#include <vector>

using namespace std;

typedef struct {
    char *key;
    char *value;
} rsArgument;

namespace rstools {
namespace batch {
namespace util {

class RSTask {
    
    public:
        
        RSTask(const char*, const char*);
        ~RSTask();
        
        void setDescription(char*);
        char* getDescription();
        void addArgument(rsArgument*);
        vector<rsArgument*> getArguments();
        rsArgument* getArgument(char* name);
        char* getOutputPath();
        void setOutputPath(char*);
        void setShowOutput(bool showOutput);
        bool shouldShowOutput();
        char** getCallString(int *argc);
        const char* getName();
        const char* getCode();
        void setCmd(char* cmd);
        char* getCmd();
        
        static RSTask* taskFactory(const char *code);
    
    protected:
        const char* name;
        const char* code;
        char *description;
        vector<rsArgument*> arguments;
        char *outputPath;
        char *cmd;
        bool showOutput;
};

}}} // namespace rstools::batch::util

#endif
