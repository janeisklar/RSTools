#include "rsjob.hpp"

namespace rstools {
namespace batch {
namespace util {

RSJob::RSJob(char* jobfile)
{
    setJobfile(jobfile);
}

RSJob::~RSJob()
{
    
}

void RSJob::setJobfile(char* jobfile)
{
    this->jobfile = jobfile;
} 

char* RSJob::getJobfile() {
    return this->jobfile;
}

void RSJob::setDescription(char* description)
{
    this->description = description;
}

char* RSJob::getDescription() {
    return this->description;
}

void RSJob::addTask(RSTask* task) {
    this->tasks.push_back(task);
}

vector<RSTask*> RSJob::getTasks() {
    return this->tasks;
}

void RSJob::addArgument(rsArgument* argument)
{
    this->arguments.push_back(argument);
}

vector<rsArgument*> RSJob::getArguments()
{
    return this->arguments;
}

char* RSJob::toXml()
{
    char *oldXml;
    char *xml = rsStringConcat(
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
        "<job xmlns=\"http://www.fmri.at/rstools\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.fmri.at/rstools http://www.fmri.at/rstools/job.xsd\">\n",
        "    <parameters>",
        NULL
    );
    
    for ( vector<rsArgument*>::size_type j = 0; j != arguments.size(); j++ ) {
        rsArgument *arg = arguments[j];
        oldXml = xml;
        char *argXml = this->_argumentToXml(arg);
        xml = rsStringConcat(xml, argXml, NULL);
        rsFree(argXml);
        rsFree(oldXml);
    }
    
    
    oldXml = xml;
    xml = rsStringConcat(
        xml,
        "\n",
        "    </parameters>\n",
        "    <tasks>\n",
        NULL
    );
    rsFree(oldXml);
    
    for ( vector<RSTask*>::size_type j = 0; j != tasks.size(); j++ ) {
        RSTask *task = tasks[j];
        char *oldXml = xml;
        char *taskXml = task->toXml();
        xml = rsStringConcat(xml, taskXml, NULL);
        rsFree(taskXml);
        rsFree(oldXml);
    }
    
    oldXml = xml;
    xml = rsStringConcat(
        xml, 
        "    </tasks>\n",
        "</job>\n",
        NULL
    );
    rsFree(oldXml);
    
    return xml;
}

char* RSJob::_argumentToXml(rsArgument *arg)
{
    return rsStringConcat(
        "\n",
        "        <param name=\"", arg->key, "\">\n",
        "            ", arg->value, "\n",
        "        </param>",
        NULL);
}
}}} // namespace rstools::batch::util
