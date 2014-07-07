#include "rsjob.hpp"

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
