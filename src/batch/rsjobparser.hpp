#ifndef rstools_rsjobparser_hpp
#define rstools_rsjobparser_hpp

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <iostream>
#include <vector>
#include "rsjob.hpp"
#include "src/rscommon.h"
#include "parseErrorHandler.hpp"

using namespace std;
using namespace xercesc;

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
		
	protected:
		void parseParameters();
		void parseTasks();
		void parseTask();
		
		char *leftTrimString(char *s);
		char *rightTrimString(char *s);
		char *trimString(char *s);
		char *replaceString(char *orig, char *rep, char *with);
};

#endif
