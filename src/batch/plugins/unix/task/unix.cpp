#include "unix.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {
namespace task {
   
Unix::Unix(const char* code, const char* name) : RSUnixTask(code, name)
{
    setCmd((char*)"");
}
    
void Unix::fillInJobArguments(RSJob* job, RSJobParser* parser)
{
    vector<rsArgument*> jobArguments = job->getArguments();
    
    for ( vector<rsArgument*>::size_type i = 0; i < jobArguments.size(); i++ ) {
        rsArgument *arg = jobArguments[i];
        
        // create placeholder string
        char *key = (char*)malloc((strlen(arg->key)+(size_t)4)*sizeof(char));
        sprintf(key, "${%s}", arg->key);
    
        // replace
        setCmd(
            parser->replaceString(getCmd(false), key, arg->value)
        );
        
        free(key);
    }
}

void Unix::parseTaskFromXml(DOMNodeIterator* walker, DOMNode* &current_node)
{
    for ( current_node = walker->nextNode(); current_node != NULL; current_node = walker->nextNode() ) {

        char *thisNodeName   = XMLString::transcode(current_node->getNodeName());
        char *parentNodeName = XMLString::transcode(current_node->getParentNode()->getNodeName());

        if( strcmp(parentNodeName, getCode()) ) {
            break;
        }
        
        if( ! strcmp(thisNodeName, "description") ) {
                char *desc = XMLString::transcode(current_node->getFirstChild()->getNodeValue());
                setDescription(desc);
        } else if( ! strcmp(thisNodeName, "cmd") ) {
            char *cmd = rsTrimString(XMLString::transcode(current_node->getFirstChild()->getNodeValue()));
            setCmd(cmd);
        }
    }
}

void Unix::setCmd(const char* cmd)
{   
    // set argument as well
    const char* optionName = "command";
    rsArgument* argument = getArgument(optionName);
    
    if ( argument == NULL ) {
          argument = (rsArgument*)malloc(sizeof(rsArgument));
          argument->key = rsString(optionName);
          addArgument(argument);
    }
    
    argument->value = rsString(cmd);
} 

char* Unix::getCmd(bool asExecuted) {
    rsArgument* argument = getArgument("command");
    
    if ( argument == NULL ) {
        return NULL;
    }
    
    return argument->value;
}

char* Unix::toXml()
{
    return rsStringConcat(
        "        <", this->code, ">\n",
        "            <description>", this->getDescription(), "</description>\n",
        "            <cmd>\n",
        this->getCmd(false),
        "\n",
        "            </cmd>\n",
        "        </", this->code, ">\n",
        NULL
    );
}
}}}}} // namespace rstools::batch::plugins::unix::task
