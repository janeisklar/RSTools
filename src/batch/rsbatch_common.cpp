#include "rsbatch_common.hpp"

#include <cstdlib>
#include <unistd.h>

using namespace std;
using namespace rstools::batch::util;

void rsBatchInit(rsBatchParameters* p)
{
    p->parametersValid = FALSE;
    p->verbose = p->quiet == FALSE;

    // ensure that plugins are loaded
    PluginManager::getInstance().loadPlugins();
    
    if ( p->jobpath == NULL ) {
        fprintf(stderr, "Error: No jobfile was specified!\n");
        return;
    }
    
    // Parse job-file
    RSJobParser* parser = new RSJobParser(p->jobpath);
    
    if ( ! parser->parse() ) {
        fprintf(stderr, "Job-file could not be read!\n");
        return;
    }
    
    // Fill in user arguments into job arguments
    p->job = parser->getJob();
    parser->fillInUserArguments(p->arguments, p->nArguments);
    
    // Check if any arguments are missing
    vector<char*> missingArguments = parser->getMissingArguments();

    if ( p->viewArgument == NULL && missingArguments.size() > 0 ) {
        fprintf(stderr, "The following arguments need to be supplied to run the job-file!\n");
        for ( vector<char*>::size_type t = 0; t < missingArguments.size(); t++ ) {
            char *argument = missingArguments[t];
            fprintf(stderr, " * %s\n", argument);
        }
        return;
    }
    
    p->parametersValid = TRUE;
}

void rsBatchRun(rsBatchParameters *p)
{   
    // Check if all that needs to be done is to print a specified parameter
    if ( p->viewArgument != NULL ) {
        return rsBatchPrintParameter(p);
    }
    
    // Parse the runtime arguments
    vector<RSTask*> tasks = p->job->getTasks();
    RSTool** tools = (RSTool**)rsMalloc(tasks.size()*sizeof(RSTool*));
    
    // Prepare tasks    
    for ( vector<RSTask*>::size_type t = 0; t < tasks.size(); t++ ) {
        
        if ( rsBatchTaskShouldBeSkipped(p, (int)t) ) {
            continue;
        }
        
        RSTask *task = tasks[t];
        
        // create the tool that's needed to execute the current task
        RSTool *tool = RSTool::toolFactory(task->getCode());
        tool->setTask(task);
        tool->setThreads(p->threads);
        
        // create the appropriate call string just like we'd be calling it from the command line
        int argc;
        char **argv = task->getCallString(&argc);

        // let the tool parse the call string
        tool->parseParams(argc, argv);
        
        // check if the tool is happy and able to continue with the arguments it received
        if ( ! tool->isEverythingFine() ) {
            return rsBatchPrintExecutionError(tool, (t+1), "parsing the supplied arguments");
        }
        
        tools[t] = tool;        
    }
    
    // Check if instead of executing only an overview over all tasks should be given
    if ( p->showOverview ) {
        return rsBatchShowJobOverview(p, tools);
    }
    
    if ( p->verbose ) {
        fprintf(stdout, "# Executing tasks\n");
    }
    
    // Execute tasks
    for ( short t = 0; t < tasks.size(); t++ ) {
        
        if ( rsBatchTaskShouldBeSkipped(p, (int)t) ) {
            continue;
        }
        
        RSTool *tool = tools[t];

        if ( p->verbose ) {
            fprintf(stdout, "Task %d of %d\n", (int)t+1, (int)tasks.size());
            fprintf(stdout, "Description: %s\n", tool->getTask()->getDescription());
            fprintf(stdout, ">Initializing..\n");
        }

        // initialize tool
        tool->init();

        // check for errors during initialization
        if ( ! tool->isEverythingFine() ) {
            return rsBatchPrintExecutionError(tool, (t+1), "initializing");
        }

        if ( p->verbose ) {
            fprintf(stdout, ">Running..\n");
        }

        // run tool
        tool->run();

        // check for errors during execution
        if ( ! tool->isEverythingFine() ) {
            return rsBatchPrintExecutionError(tool, (t+1), "running");
        }
        
        // cleanup
        tool->destroy();

        if ( p->verbose ) {
            fprintf(stdout, "\n\n");
        }
    }
}

void rsBatchPrintExecutionError(RSTool *tool, int taskNum, char const * state)
{
    fprintf(stderr, "\nError in task #%d while %s\n", taskNum, state);
    tool->printCallString(stderr);
    
    if ( strcmp(tool->getOutput(), "") ) {
        fprintf(stderr, "\n\nMessage:\n%s\n", tool->getOutput());
    }
}

void rsBatchPrintParameter(rsBatchParameters *p)
{
    vector<rsArgument*> arguments = p->job->getArguments();
    for ( vector<rsArgument*>::size_type i = 0; i != arguments.size(); i++ ) {
        rsArgument *arg = arguments[i];
        if ( ! strcmp(p->viewArgument, arg->key) ) {
            fprintf(stdout, "%s", arg->value);
            return;
        }
    }
    
    fprintf(stdout, "Parameter '%s' is not contained within jobfile '%s'!\n", p->viewArgument, p->jobpath);
}

void rsBatchDestroy(rsBatchParameters* p)
{
    // Destruct the XML Parser
    XMLPlatformUtils::Terminate();
}

void rsBatchShowJobOverview(rsBatchParameters* p, RSTool** tools) {
    
    vector<RSTask*>::size_type nTasks = p->job->getTasks().size();
    
    for ( vector<RSTask*>::size_type t=0; t<nTasks; t++ ) {
        
        if ( rsBatchTaskShouldBeSkipped(p, (int)t) ) {
            continue;
        }
        
        RSTool *tool = tools[t];
        RSTask *task = tool->getTask();
        
        char **chunkedDescription;
        size_t nLines;
        
        rsStringWordWrap(task->getDescription(), &chunkedDescription, &nLines, 57);
        
        fprintf(stdout, "################################################################\n");
        fprintf(stdout, "## Task #%-3d                                                  ##\n", (int)(t+1));
        for ( size_t i=0; i<nLines; i++ ) {
            fprintf(stdout, "## %-58s ##\n", chunkedDescription[i]);
        }
        fprintf(stdout, "################################################################\n");
        fprintf(stdout, "\n");
        
        tool->printCallString(stdout);
    }
}

BOOL rsBatchTaskShouldBeSkipped(rsBatchParameters* p, int taskId)
{
    size_t t = 0;
    while ( p->tasksToSkip[t] != -1 ) {
        
        if ( p->tasksToSkip[t] == taskId ) {
            return TRUE;
        }
        
        t++;
    }
    
    return FALSE;
}
