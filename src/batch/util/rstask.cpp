#include "rstask.hpp"
#include "rstool.hpp"

namespace rstools {
namespace batch {
namespace util {

RSTask::RSTask(const char* code, const char* name)
{
    this->code = code;
    this->name = name;
    this->job = NULL;
    setOutputPath(NULL);
    setDescription((char*)"");
    setShowOutput(false);
}

RSTask::~RSTask()
{}

void RSTask::setDescription(char* description)
{
    this->description = description;
}

char* RSTask::getDescription()
{
    return this->description;
}

void RSTask::setShowOutput(bool showOutput)
{
    this->showOutput = showOutput;
}

bool RSTask::shouldShowOutput()
{
    return this->showOutput;
}

RSJob* RSTask::getJob()
{
    return this->job;
}

void RSTask::addArgument(rsArgument* argument)
{
    this->arguments.push_back(argument);
}

vector<rsArgument*> RSTask::getArguments()
{
    return this->arguments;
}

rsArgument* RSTask::getArgument(const char* name)
{
    vector<rsArgument*> arguments = this->getArguments();
    for ( vector<rsArgument*>::size_type i = 0; i != arguments.size(); i++ ) {
        rsArgument *arg = arguments[i];
        
        if ( ! strcmp(arg->key, name) ) {
            return  arg;
        }
    }
    
    return NULL;
}

/*
 * Returns the value of an argument if it was supplied or a 
 * default value from the config instead
 */
const char* RSTask::getDefaultArgumentValue(const char* name)
{
    rsArgument* arg = getArgument(name);
    if ( arg != NULL ) {
        return arg->value;
    }
    
    arg = RSConfig::getInstance().getArgument(name);
    
    if ( arg != NULL ) {
        return arg->value;
    }
    
    return NULL;
}

void RSTask::removeArgument(const char* name)
{
    vector<rsArgument*> arguments = this->getArguments();
    for ( vector<rsArgument*>::size_type i = 0; i != arguments.size(); i++ ) {
        rsArgument *arg = arguments[i];
        
        if ( ! strcmp(arg->key, name) ) {
            arguments.erase(arguments.begin() + i);
            return;
        }
    }
}

char* RSTask::getOutputPath()
{
    return this->outputPath;
}

void RSTask::setOutputPath(char* outputPath)
{
    this->outputPath = outputPath;
}

RSTask* RSTask::taskFactory(const char *code)
{
    rsToolRegistration* registration = RSTool::findRegistration(code);
    
    if ( registration == NULL ) {
        throw std::invalid_argument(string("Task with the code '") + string(code) + string("' unknown"));
    }
    
    return registration->createTask();
}

/*
 * Returns the string it would be called with when run outside of the pipleline
 * from within bash
 */
char** RSTask::getCallString(int *argc)
{
    vector<rsArgument*> arguments = this->getArguments();
    char **argv = (char**)malloc(sizeof(char*)*(arguments.size()+1));
    *argc = (int)arguments.size()+1;
    
    const char* taskName = this->getCode();
    argv[0] = (char*)malloc((strlen(taskName)+1)*sizeof(char));
    sprintf(argv[0], "%s", taskName);
    
    for ( vector<rsArgument*>::size_type i = 0; i != arguments.size(); i++ ) {
        rsArgument *arg = arguments[i];
        
        if ( arg->value == NULL ) { // --key
            size_t length = strlen(arg->key) + 3;
            argv[i+1] = (char*)malloc(length*sizeof(char));
            sprintf(argv[i+1], "--%s", arg->key);
        } else { // --key=value
            size_t length = strlen(arg->key) + strlen(arg->value) + 4;
            argv[i+1] = (char*)malloc(length*sizeof(char));
            sprintf(argv[i+1], "--%s=%s", arg->key, arg->value);
        }
    }
    
    return argv;
}

const char* RSTask::getName()
{
    return name;
}

const char* RSTask::getCode()
{
    return code;
}

void RSTask::fillInJobArguments(RSJob* job, RSJobParser* parser)
{
    this->job = job;
    vector<rsArgument*> jobArguments = job->getArguments();
    
    for ( vector<rsArgument*>::size_type i = 0; i < jobArguments.size(); i++ ) {
        rsArgument *arg = jobArguments[i];

        for ( vector<rsArgument*>::size_type j = 0; j != arguments.size(); j++ ) {

            // create placeholder string
            char *key = (char*)malloc((strlen(arg->key)+(size_t)4)*sizeof(char));
            sprintf(key, "${%s}", arg->key);

            // replace
            rsArgument *arg2 = arguments[j];                
            arg2->value = parser->replaceString(arg2->value, key, arg->value);

            free(key);
        }
    }
}

void RSTask::parseTaskFromXml(DOMNodeIterator* walker, DOMNode* &current_node)
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
        } else if( ! strcmp(thisNodeName, "args") ) {
            // parse arguments
            for ( current_node = walker->nextNode(); current_node != NULL; current_node = walker->nextNode() ) {
                
                parentNodeName = XMLString::transcode(current_node->getParentNode()->getNodeName());

                if( strcmp(parentNodeName, "args") ) {
                    break;
                }
                
                rsArgument *arg = (rsArgument*)malloc(sizeof(rsArgument));
                arg->key   = rsTrimString(XMLString::transcode(current_node->getAttributes()->getNamedItem(XMLString::transcode("name"))->getNodeValue()));
                DOMNode *valueNode = current_node->getFirstChild();
                if ( valueNode == NULL ) {
                    arg->value = NULL;
                } else {
                    arg->value = rsTrimString(XMLString::transcode(valueNode->getNodeValue()));
                }
                addArgument(arg);
            }
            current_node = walker->previousNode();
        } else if( ! strcmp(thisNodeName, "options") ) {
            // parse options
            for ( current_node = walker->nextNode(); current_node != NULL; current_node = walker->nextNode() ) {
                
                thisNodeName = XMLString::transcode(current_node->getNodeName());
                parentNodeName = XMLString::transcode(current_node->getParentNode()->getNodeName());
                
                if( strcmp(parentNodeName, "options") ) {
                    break;
                }
                
                char *value = rsTrimString(XMLString::transcode(current_node->getFirstChild()->getNodeValue()));
                                
                if ( ! strcmp(thisNodeName, "save_output") ) {
                    setOutputPath(value);
                } else if ( ! strcmp(thisNodeName, "show_output") ) {
                    setShowOutput( ! strcmp(value, "1") );
                }
            }
            current_node = walker->previousNode();
        }
    }
}

char* RSTask::toXml()
{
    char *xml = rsStringConcat(
        "        <", this->code, ">\n",
        "            <description>", this->getDescription(), "</description>\n",
        "            <args>",
        NULL
    );
    
    for ( vector<rsArgument*>::size_type j = 0; j != arguments.size(); j++ ) {
        rsArgument *arg = arguments[j];
        
        if ( arg->value != NULL && !strcmp(arg->value,"") ) {
            continue;
        }
        
        char *oldXml = xml;
        char *argXml = this->_argumentToXml(arg);
        xml = rsStringConcat(xml, argXml, NULL);
        rsFree(argXml);
        rsFree(oldXml);
    }
    
    char *oldXml = xml;
    xml = rsStringConcat(
        xml, 
        "\n",
        "            </args>\n",
        "        </", this->code, ">\n",
        NULL
    );
    rsFree(oldXml);
    
    return xml;
}

char* RSTask::_argumentToXml(rsArgument *arg)
{
    if ( arg->value == NULL ) {
        return rsStringConcat("\n                <arg name=\"", arg->key, "\"/>", NULL);
    } else {
        return rsStringConcat("\n                <arg name=\"", arg->key, "\">", arg->value, "</arg>", NULL);
    }
}

}}} // namespace rstools::batch::util
