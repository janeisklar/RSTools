#ifndef rstools_rsbatch_util_rsjobparser_hpp
#define rstools_rsbatch_util_rsjobparser_hpp

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include "rstool.hpp"
#include "rsjob.hpp"
#include "src/rscommon.h"
#include "src/utils/rsstring.h"
#include "parseErrorHandler.hpp"

using namespace std;
using namespace xercesc;

namespace rstools {
namespace batch {
namespace util {

class RSJob;
    
class RSJobParser {
    protected:
        RSJob*           job;
        XercesDOMParser* parser;
        ErrorHandler*    errHandler;
        DOMDocument*     doc;
        DOMElement*      docRootNode;
        DOMNodeIterator* walker;
        DOMNode*         current_node;
        
    public:
        RSJobParser(char* jobfile);
        ~RSJobParser();
        bool parse();
        RSJob* getJob();
        void fillInUserArguments(rsArgument **arguments, const short nArguments);
        vector<char*> getMissingArguments();
        char *replaceString(const char *orig, const char *rep, const char *with);
        
    protected:
        void parseParameters();
        void parseTasks();
        void parseTask();
        
        char *mergePluginXSDExtensions();
        
        static inline void replaceAll(std::string &str, const std::string& from, const std::string& to)
        {
            size_t start_pos = 0;
            while((start_pos = str.find(from, start_pos)) != std::string::npos) {
                str.replace(start_pos, from.length(), to);
                start_pos += to.length();
            }
        }
};

}}} // namespace rstools::batch::util

#endif
