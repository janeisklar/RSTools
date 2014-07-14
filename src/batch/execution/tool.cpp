#include "tool.hpp"

#include "timecourse.hpp"
#include "regression.hpp"
#include "bandpass.hpp"
#include "motionscrubbing.hpp"
#include "correlation.hpp"
#include "roi.hpp"
#include "unix.hpp"

namespace rstools {
namespace batch {
namespace execution {

const short Tool::tools[] = {
    RSTask::TASK_RSTIMECOURSE,
    RSTask::TASK_RSREGRESSION,
    RSTask::TASK_RSBANDPASS,
    RSTask::TASK_RSMOTIONSCRUBBING,
    RSTask::TASK_RSCORRELATION,
    RSTask::TASK_RSROI,
    RSTask::TASK_UNIX,
    -1
};
    
Tool* Tool::factory(const short taskCode) {

    switch ( taskCode ) {
        case RSTask::TASK_RSTIMECOURSE:
            return new Timecourse();
        case RSTask::TASK_RSREGRESSION:
            return new Regression();
        case RSTask::TASK_RSBANDPASS:
            return new Bandpass();
        case RSTask::TASK_RSMOTIONSCRUBBING:
            return new MotionScrubbing();
        case RSTask::TASK_RSCORRELATION:
            return new Correlation();
        case RSTask::TASK_RSROI:
            return new Roi();
        case RSTask::TASK_UNIX:
            return new Unix();
    }
    
    throw std::invalid_argument("taskCode unknown");
}

void Tool::parseParams(int argc, char** argv)
{
    oc = new rstools::batch::OutputCatcher();
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

void Tool::init()
{
    oc = new rstools::batch::OutputCatcher();
    oc->beginCapture();

    this->_init();
    
    oc->endCapture();
    _message = oc->getCapture();
    delete oc;
}

void Tool::run()
{
    if ( this->getTask()->shouldShowOutput() ) {
        this->_run();
    } else {
        oc = new rstools::batch::OutputCatcher();
        oc->beginCapture();
    
        this->_run();

        oc->endCapture();
        _message = oc->getCapture();
        delete oc;
    }
}

char const* Tool::getOutput()
{
    return _message.c_str();
}

char** Tool::getCallString(int *argc)
{
    *argc = this->argc;
    return this->argv;
}

RSTask* Tool::getTask()
{
    return this->task;
}

void Tool::setTask(RSTask *task)
{
    this->task = task;
}

void Tool::setThreads(int threads)
{
    this->threads = threads;
}

void Tool::showProgressCallback(rsReportProgressEvent *event, void *userdata)
{
    rstools::batch::OutputCatcher *oc = (rstools::batch::OutputCatcher*) userdata;
    printProgressBar(oc->getStdout(), event->percentage, event->run, (char*)"");
}

void Tool::printProgressBar(FILE* stream, double percentage, int run, char* description)
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

}}} // namespace rstools::batch::execution
