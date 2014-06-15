#include "rstask.hpp"

const short RSTask::TASK_UNIX;
const short RSTask::TASK_RSTIMECOURSE;
const short RSTask::TASK_RSREGRESSION;
const short RSTask::TASK_RSBANDPASS;
const short RSTask::TASK_RSMOTIONSCRUBBING;
const short RSTask::TASK_RSCORRELATION;
const short RSTask::TASK_RSROI;

RSTask::RSTask(const short taskCode)
{
	task = taskCode;
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

short RSTask::getTask()
{
	return this->task;
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

short RSTask::getTaskFromName(char *name)
{
	if ( !strcicmp(name, "rstimecourse") ) {
		return RSTask::TASK_RSTIMECOURSE;
	}
	
	if ( !strcicmp(name, "rsregression") ) {
		return RSTask::TASK_RSREGRESSION;
	}
	
	if ( !strcicmp(name, "rsbandpass") ) {
		return RSTask::TASK_RSBANDPASS;
	}
	
	if ( !strcicmp(name, "rsmotionscrubbing") ) {
		return RSTask::TASK_RSMOTIONSCRUBBING;
	}
	
	if ( !strcicmp(name, "rscorrelation") ) {
		return RSTask::TASK_RSCORRELATION;
	}
	
	if ( !strcicmp(name, "rsroi") ) {
		return RSTask::TASK_RSROI;
	}
	
	if ( !strcicmp(name, "unix") ) {
		return RSTask::TASK_UNIX;
	}
	
	return RSTask::TASK_UNKNOWN;
}

char const* RSTask::getNameForTask(short code)
{
	switch ( code ) {
		case RSTask::TASK_RSTIMECOURSE:
			return "rstimecourse";
		case RSTask::TASK_RSREGRESSION:
			return "rsregression";
		case RSTask::TASK_RSBANDPASS:
			return "rsbandpass";
		case RSTask::TASK_RSMOTIONSCRUBBING:
			return "rsmotionscrubbing";
		case RSTask::TASK_RSCORRELATION:
			return "rscorrelation";
		case RSTask::TASK_RSROI:
			return "rsroi";
		case RSTask::TASK_UNIX:
			return "unix";
		default:
			return "unknown";
	}
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

char const * RSTask::getName()
{
	return RSTask::getNameForTask(this->getTask());
}

int RSTask::strcicmp(char const *a, char const *b)
{
    for (;; a++, b++) {
        int d = tolower(*a) - tolower(*b);
        if (d != 0 || !*a)
            return d;
    }
}