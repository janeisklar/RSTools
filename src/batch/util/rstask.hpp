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
        
        virtual void setDescription(char*);
        virtual char* getDescription();
        virtual void addArgument(rsArgument*);
        virtual vector<rsArgument*> getArguments();
        virtual rsArgument* getArgument(const char* name);
        virtual void removeArgument(const char* name);
        virtual RSJob* getJob();
        virtual char* getOutputPath();
        virtual void setOutputPath(char*);
        virtual void setShowOutput(bool showOutput);
        virtual bool shouldShowOutput();
        virtual char** getCallString(int *argc);
        virtual const char* getName();
        virtual const char* getCode();
        virtual void fillInJobArguments(RSJob* job, RSJobParser* parser);
        virtual void parseTaskFromXml(DOMNodeIterator* walker, DOMNode* &current_node);
        virtual char* toXml();
        virtual char* _argumentToXml(rsArgument *arg);
        
        static RSTask* taskFactory(const char *code);
    
    protected:
        const char* name;
        const char* code;
        char *description;
        vector<rsArgument*> arguments;
        char *outputPath;
        bool showOutput;
        RSJob *job;
};

}}} // namespace rstools::batch::util

#endif
