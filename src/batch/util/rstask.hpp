#ifndef rstools_rsbatch_util_rstask_hpp
#define rstools_rsbatch_util_rstask_hpp

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char *key;
    char *value;
} rsArgument;

#ifdef __cplusplus
}
#endif

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <locale>
#include <vector>
#include "rsjob.hpp"
#include "rsjobparser.hpp"

using namespace std;
using namespace xercesc;

namespace rstools {
namespace batch {
namespace util {

class RSJob;
class RSJobParser;

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
        void fillInJobArguments(RSJob* job, RSJobParser* parser);
        void parseTaskFromXml(DOMNodeIterator* walker, DOMNode* &current_node);
        
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
