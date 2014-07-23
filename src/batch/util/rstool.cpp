#include "rstool.hpp"

namespace rstools {
namespace batch {
namespace util {

vector<rsToolRegistration*> RSTool::tools;
vector<rsXSDExtension*> RSTool::xsdExtensions;

rsToolRegistration* RSTool::findRegistration(const char* code)
{
    for(vector<rsToolRegistration*>::iterator it = tools.begin(); it != tools.end(); ++it) {
        rsToolRegistration* registration = *it;
        
        if ( ! strcmp(registration->code, code) ) {
            return registration;
        }
    }
    
    return NULL;
}
    
RSTool* RSTool::toolFactory(const char* code)
{
    rsToolRegistration *registration = findRegistration(code);
    
    if ( registration == NULL ) {
        throw std::invalid_argument(string("tool '") + string(code) + string("' is unknown"));
    }
    
    return registration->createTool();
}

void RSTool::registerTool(rsToolRegistration* registration)
{
    tools.push_back(registration);
}


void RSTool::registerXSDExtension(rsXSDExtension* extension)
{
    xsdExtensions.push_back(extension);
}

vector<rsXSDExtension*> RSTool::getXSDExtensions()
{
    return xsdExtensions;
}

vector<const char*> RSTool::getTools()
{
    vector<const char*> toolCodes;
    for(vector<rsToolRegistration*>::iterator it = tools.begin(); it != tools.end(); ++it) {
        rsToolRegistration* registration = *it;
        toolCodes.push_back(registration->code);
    }
    return toolCodes;
}

void RSTool::parseParams(int argc, char** argv)
{
    oc = new OutputCatcher();
    oc->beginCapture();
    
    // copy call arguments
    this->argc = argc;
    this->argv = (char**)malloc(argc*sizeof(char*));
    
    for ( int i=0; i<argc; i++ ) {
        this->argv[i] = (char*)malloc((strlen(argv[i])+1)*sizeof(char));
        sprintf(this->argv[i], "%s", argv[i]);
    }
    
    // parse them
    this->_parseParams(argc, argv);
    
    oc->endCapture();
    _message = oc->getCapture();
    delete oc;
}

void RSTool::init()
{
    oc = new OutputCatcher();
    oc->beginCapture();

    this->_init();
    
    oc->endCapture();
    _message = oc->getCapture();
    delete oc;
}

void RSTool::run()
{
    if ( this->getTask()->shouldShowOutput() ) {
        this->_run();
    } else {
        oc = new OutputCatcher();
        oc->beginCapture();
    
        this->_run();

        oc->endCapture();
        _message = oc->getCapture();
        delete oc;
    }
}

char const* RSTool::getOutput()
{
    return _message.c_str();
}

char** RSTool::getCallString(int *argc)
{
    *argc = this->argc;
    return this->argv;
}

RSTask* RSTool::getTask()
{
    return this->task;
}

void RSTool::setTask(RSTask *task)
{
    this->task = task;
}

void RSTool::setThreads(int threads)
{
    this->threads = threads;
}

void RSTool::showProgressCallback(rsReportProgressEvent *event, void *userdata)
{
    OutputCatcher *oc = (OutputCatcher*) userdata;
    printProgressBar(oc->getStdout(), event->percentage, event->run, (char*)"");
}

void RSTool::printProgressBar(FILE* stream, double percentage, int run, char* description)
{
    // determine the width of the terminal window
    struct winsize max;
    ioctl(0, TIOCGWINSZ , &max);

    const int width = max.ws_col;
    const int progressBarWidth = width - 9;

    // remove old progress bar if present
    if ( run > 0 ) {
        for ( int i=0; i<5*width; i++ ) {
            fprintf(stream, "\b");
        }
    }

    // print description
    fprintf(stream, "%s", description);

    // fake new line
    for ( int i = strlen(description); (i % width) != 0; i++ ) {
        fprintf(stream, " ");
    }

    fprintf(stream, "[");
    const int done = (int)ceil((percentage / 100.0) * (double)progressBarWidth);
    const int missing = progressBarWidth - done;

    for ( int i=0; i<done; i++ ) {
        fprintf(stream, "+");
    }

    for ( int i=0; i<missing; i++ ) {
        fprintf(stream, "-");
    }

    fprintf(stream, "] % 3d%%", (int)ceil(percentage));
}

void RSTool::printCallString(FILE *stream)
{
    int argc;
    char **argv = getCallString(&argc);

    fprintf(stream, "Tool:\n %s\n\nParams:\n", getTask()->getName());
    for ( int i=1; i<argc; i++ ) {
        fprintf(stream, "  %s\n", argv[i]);
    }
    
    fprintf(stream, "\n");
}

}}} // namespace rstools::batch::util
