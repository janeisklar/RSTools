#ifndef rstools_parseerrorhandler_hpp
#define rstools_parseerrorhandler_hpp

#include <stdio.h>
 
#include <xercesc/util/XMLString.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/validators/common/Grammar.hpp>

XERCES_CPP_NAMESPACE_USE

class ParserErrorHandler : public ErrorHandler
{
    private:
        void reportParseException(const SAXParseException& ex)
        {
            char* msg = XMLString::transcode(ex.getMessage());
			
            fprintf(stderr, "Parsing problem at column %lu, line %lu, %s\n", (long unsigned int)ex.getColumnNumber(), (long unsigned int)ex.getLineNumber(), msg);
            XMLString::release(&msg);
			throw ex;
        }
 
    public:
        void warning(const SAXParseException& ex)
        {
            reportParseException(ex);
        }
 
        void error(const SAXParseException& ex)
        {
            reportParseException(ex);
        }
 
        void fatalError(const SAXParseException& ex)
        {
            reportParseException(ex);
        }
 
        void resetErrors()
        {
        }
};

#endif