#ifndef rstools_rsbatch_execution_tool_h
#define rstools_rsbatch_execution_tool_h

#include <iostream>
#include "src/batch/outputCatcher.hpp"
#include "src/batch/rstask.hpp"
#include <sys/ioctl.h>
#include <math.h>
#include "src/rscommon.h"

namespace rstools {
namespace batch {
namespace execution {

class Tool {

    public:
        
        void parseParams(int argc, char** argv)
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
        
        void init()
        {
            oc = new rstools::batch::OutputCatcher();
            oc->beginCapture();

            this->_init();
            
            oc->endCapture();
            _message = oc->getCapture();
            delete oc;
        }
        
        void run()
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
        
        char const * getOutput()
        {
            return _message.c_str();
        }
        
        char **getCallString(int *argc)
        {
            *argc = this->argc;
            return this->argv;
        }
        
        RSTask *getTask()
        {
            return this->task;
        }
        
        void setTask(RSTask *task)
        {
            this->task = task;
        }
        
        void setThreads(int threads)
        {
            this->threads = threads;
        }
        
        static void showProgressCallback(rsReportProgressEvent *event, void *userdata)
        {
            rstools::batch::OutputCatcher *oc = (rstools::batch::OutputCatcher*) userdata;
            printProgressBar(oc->getStdout(), event->percentage, event->run, (char*)"");
        }
        
        virtual void destroy() = 0;
        virtual bool isEverythingFine() = 0;
        
        static void printProgressBar(FILE* stream, double percentage, int run, char* description)
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
        
    protected:
        virtual void _parseParams(int, char**) = 0;
        virtual void _init() = 0;
        virtual void _run() = 0;
        std::string _message;
        int argc;
        char** argv;
        RSTask *task;
        int threads;
        rstools::batch::OutputCatcher *oc;
};

}}} // namespace rstools::batch::execution

#endif
