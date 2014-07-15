#include "rstask.hpp"
#include "rstool.hpp"

namespace rstools {
namespace batch {
namespace util {

RSTask::RSTask(const char* code, const char* name)
{
    this->code = code;
    this->name = name;
    setOutputPath(NULL);
    setDescription((char*)"");
    setCmd((char*)"");
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

void RSTask::setCmd(char* cmd)
{
    this->cmd = cmd;
} 

char* RSTask::getCmd() {
    return this->cmd;
}

void RSTask::addArgument(rsArgument* argument)
{
    this->arguments.push_back(argument);
}

vector<rsArgument*> RSTask::getArguments()
{
    return this->arguments;
}

char* RSTask::getOutputPath()
{
    return this->outputPath;
}

void RSTask::setOutputPath(char* outputPath)
{
    this->outputPath = outputPath;
}

RSTask* RSTask::taskFactory(const char *name)
{
    rsToolRegistration* registration = RSTool::findRegistration(name);
    
    if ( registration == NULL ) {
        throw std::invalid_argument(string("Task '") + string(name) + string("' unknown"));
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
    
    const char* taskName = this->getName();
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

}}} // namespace rstools::batch::util
