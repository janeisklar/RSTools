#include "rsjobparser.hpp"

RSJobParser::RSJobParser(char *jobfile)
{
	this->job = new RSJob(jobfile);
	this->errHandler = (ErrorHandler*) new ParserErrorHandler();
}

RSJobParser::~RSJobParser()
{
	delete this->job;
	delete this->errHandler;
}

bool RSJobParser::parse()
{
	// Initialize the XML Parser
	try {
    	XMLPlatformUtils::Initialize();
  	} catch (const XMLException& e) {
		char *message = XMLString::transcode(e.getMessage());
		fprintf(stderr, "Job File Error: %s\n", message);
		XMLString::release(&message);
    	return false;
  	}
	
	const char schemaLocationTemplate[] = "http://www.fmri.at/rstools %s/%s/jobs/job.xsd";
	char schemaLocation[strlen(schemaLocationTemplate)+strlen(RSTOOLS_DATA_DIR)+strlen(PACKAGE)];
	sprintf(&schemaLocation[0], schemaLocationTemplate, RSTOOLS_DATA_DIR, PACKAGE);
	
	parser = new XercesDOMParser();
	parser->setValidationScheme(XercesDOMParser::Val_Always);
	parser->setExternalSchemaLocation(schemaLocation);
	parser->setDoNamespaces(true);
    parser->setDoSchema(true);
    parser->setValidationConstraintFatal(true);
	parser->setErrorHandler(errHandler);
	
	// Parse it
	try {
	    parser->parse(XMLString::transcode(this->job->getJobfile()));
	} catch (const XMLException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		fprintf(stderr, "Job file parsing error #%d: %s\n", toCatch.getCode(), message);
		XMLString::release(&message);
		return false;
	} catch (const DOMException& toCatch) {
		char* message = XMLString::transcode(toCatch.msg);
		fprintf(stderr, "Job file parsing error #%d: %s\n", toCatch.code, message);
		XMLString::release(&message);
		return false;
	} catch (const exception& e) {
   	    fprintf(stderr, "Unexpected job file parsing error: %s\n", e.what());
	    return false;
	} catch (...) {
		fprintf(stderr, "Unexpected job file parsing error.\n");
		return false;
	}

	doc = parser->getDocument();
	docRootNode = doc->getDocumentElement();

	// Create the node iterator, that will walk through each element.
	try { 
		walker = doc->createNodeIterator(docRootNode, DOMNodeFilter::SHOW_ELEMENT, NULL, true); 
	} catch (const XMLException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		fprintf(stderr, "Job file parsing error(semantic) #%d: %s\n", toCatch.getCode(), message);
		XMLString::release(&message);
		return false;
	} catch (const DOMException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		fprintf(stderr, "Job file parsing error(semantic) #%d: %s\n", toCatch.code, message);
		XMLString::release(&message);
		return false;
	} catch (const exception& e) {
   	    fprintf(stderr, "Unexpected job file parsing(semantic) error: %s\n", e.what());
	    return false;
	} catch (...) {
		fprintf(stderr, "Unexpected job file parsing error.\n");
		return false;
	}

	// Some declarations
	current_node = NULL;
	char* thisNodeName;
	char* parentNodeName;

	// Iterate through all element nodes
	for ( current_node = walker->nextNode(); current_node != NULL; current_node = walker->nextNode() ) {

		thisNodeName = XMLString::transcode(current_node->getNodeName());
		parentNodeName = XMLString::transcode(current_node->getParentNode()->getNodeName());

		if( ! strcmp(parentNodeName, "job") ) {
			if( ! strcmp(thisNodeName, "description") ) {
				char *desc = XMLString::transcode(current_node->getFirstChild()->getNodeValue());
				this->job->setDescription(desc);
			}
		} else if( ! strcmp(parentNodeName, "parameters") ) {
			this->parseParameters();
			current_node = walker->previousNode();
		} else if ( ! strcmp(parentNodeName, "tasks") ) {
			this->parseTasks();
			current_node = walker->previousNode();
		}
	}
		
	return true;
}

RSJob* RSJobParser::getJob()
{
	return this->job;
}
	
void RSJobParser::parseParameters()
{
	char *parentNodeName = NULL;

	do {
		rsArgument *arg = (rsArgument*)malloc(sizeof(rsArgument));
		arg->key   = this->trimString(XMLString::transcode(current_node->getAttributes()->getNamedItem(XMLString::transcode("name"))->getNodeValue()));
		arg->value = this->trimString(XMLString::transcode(current_node->getFirstChild()->getNodeValue()));
		this->job->addArgument(arg);
		
		current_node = walker->nextNode();
		parentNodeName = NULL;
		if ( current_node != NULL ) {
			parentNodeName = XMLString::transcode(current_node->getParentNode()->getNodeName());
		}
	} while ( ! strcmp(parentNodeName, "parameters") );
}

void RSJobParser::parseTasks()
{
	do {
		this->parseTask();
	} while ( ! strcmp("tasks", XMLString::transcode(current_node->getParentNode()->getNodeName())) );
}

void RSJobParser::parseTask()
{
	char *taskName = XMLString::transcode(current_node->getNodeName());
	short taskCode = RSTask::getTaskFromName(taskName);
	RSTask *task = new RSTask(taskCode);
	
	for ( current_node = walker->nextNode(); current_node != NULL; current_node = walker->nextNode() ) {

		char *thisNodeName   = XMLString::transcode(current_node->getNodeName());
		char *parentNodeName = XMLString::transcode(current_node->getParentNode()->getNodeName());

		if( strcmp(parentNodeName, taskName) ) {
			break;
		}
		
		if( ! strcmp(thisNodeName, "description") ) {
				char *desc = XMLString::transcode(current_node->getFirstChild()->getNodeValue());
				task->setDescription(desc);
		} else if( ! strcmp(thisNodeName, "args") ) {
			// parse arguments
			for ( current_node = walker->nextNode(); current_node != NULL; current_node = walker->nextNode() ) {
				
				parentNodeName = XMLString::transcode(current_node->getParentNode()->getNodeName());

				if( strcmp(parentNodeName, "args") ) {
					break;
				}
				
				rsArgument *arg = (rsArgument*)malloc(sizeof(rsArgument));
				arg->key   = this->trimString(XMLString::transcode(current_node->getAttributes()->getNamedItem(XMLString::transcode("name"))->getNodeValue()));
				DOMNode *valueNode = current_node->getFirstChild();
				if ( valueNode == NULL ) {
					arg->value = NULL;
				} else {
					arg->value = this->trimString(XMLString::transcode(valueNode->getNodeValue()));
				}
				task->addArgument(arg);
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
				
				char *value = this->trimString(XMLString::transcode(current_node->getFirstChild()->getNodeValue()));
								
				if ( ! strcmp(thisNodeName, "save_output") ) {
					task->setOutputPath(value);
				} else if ( ! strcmp(thisNodeName, "show_output") ) {
					task->setShowOutput( ! strcmp(value, "1") );
				}
			}
			current_node = walker->previousNode();
		} else if( ! strcmp(thisNodeName, "cmd") ) {
			char *cmd = XMLString::transcode(current_node->getFirstChild()->getNodeValue());
			task->setCmd(cmd);
		}
	}
	
	this->job->addTask(task);
}

void RSJobParser::fillInUserArguments(rsArgument **userArguments, const short nUserArguments)
{
	vector<rsArgument*> jobArguments = job->getArguments();

	// merge job arguments and user arguments
	for ( int u=0; u<nUserArguments; u++ ) {
		rsArgument* userArgument = userArguments[u];
		bool found = false;

		// check if user argument already exists in the list of job arguments
		for ( vector<rsArgument*>::size_type i = 0; i < jobArguments.size(); i++ ) {
			rsArgument *arg = jobArguments[i];

			if ( strcmp(arg->key, userArgument->key) ) {
				continue;
			}

			// if it already exists use the value supplied by the user
			arg->value = userArgument->value;
			found = true;
			break;
		}

		// if it doesn't exist yet add it to the job arguments
		if ( ! found ) {
			rsArgument* newArgument = (rsArgument*)malloc(sizeof(rsArgument));
			newArgument->key   = (char*)malloc((strlen(userArgument->key)  +(size_t)1)*sizeof(char));
			newArgument->value = (char*)malloc((strlen(userArgument->value)+(size_t)1)*sizeof(char));
			sprintf(newArgument->key,   "%s", userArgument->key);
			sprintf(newArgument->value, "%s", userArgument->value);
			this->job->addArgument(newArgument);
		}
	}

	jobArguments = job->getArguments();

	// replace all placeholders within the job arguments
	for ( vector<rsArgument*>::size_type i = 0; i < jobArguments.size(); i++ ) {
		for ( vector<rsArgument*>::size_type j = 0; j < jobArguments.size(); j++ ) {

			// skip if both arguments are the same
			if ( i == j ) {
				continue;
			}

			rsArgument *arg = jobArguments[i];
			rsArgument *arg2 = jobArguments[j];

			// create placeholder string
			char *key = (char*)malloc((strlen(arg->key)+(size_t)4)*sizeof(char));
			sprintf(key, "${%s}", arg->key);

			// replace
			arg2->value = this->replaceString(arg2->value, key, arg->value);

			free(key);
		}
	}

	// replace all job argument placeholders within the task arguments
	vector<RSTask*> tasks = this->job->getTasks();
	for ( vector<RSTask*>::size_type t = 0; t < tasks.size(); t++ ) {
		RSTask *task = tasks[t];
		vector<rsArgument*> arguments = task->getArguments();

		for ( vector<rsArgument*>::size_type i = 0; i < jobArguments.size(); i++ ) {
			rsArgument *arg = jobArguments[i];

			for ( vector<rsArgument*>::size_type j = 0; j != arguments.size(); j++ ) {

				// create placeholder string
				char *key = (char*)malloc((strlen(arg->key)+(size_t)4)*sizeof(char));
				sprintf(key, "${%s}", arg->key);

				// replace
				rsArgument *arg2 = arguments[j];				
				arg2->value = this->replaceString(arg2->value, key, arg->value);

				free(key);
			}
		}
	}
	
	// finally, replace all job argument placeholders within the unix cmd
	
	
	for ( vector<rsArgument*>::size_type i = 0; i < jobArguments.size(); i++ ) {
		rsArgument *arg = jobArguments[i];
		
		// create placeholder string
		char *key = (char*)malloc((strlen(arg->key)+(size_t)4)*sizeof(char));
		sprintf(key, "${%s}", arg->key);
	
		for ( vector<RSTask*>::size_type t = 0; t < tasks.size(); t++ ) {
			RSTask *task = tasks[t];
			
			// replace
			task->setCmd(
				this->replaceString(task->getCmd(), key, arg->value)
			);
		}
		
		free(key);
	}
}

vector<char*> RSJobParser::getMissingArguments()
{
	vector<char*> allMissingArguments;
	
	// run throught all arguments of all tasks
	vector<RSTask*> tasks = this->job->getTasks();
	for ( vector<RSTask*>::size_type t = 0; t < tasks.size(); t++ ) {
		RSTask *task = tasks[t];
		vector<rsArgument*> arguments = task->getArguments();

		for ( vector<rsArgument*>::size_type i = 0; i < arguments.size(); i++ ) {
			char *value = arguments[i]->value;
			
			if ( value == NULL ) {
				continue;
			}
			
			size_t length = strlen(value);
			
			int start = -1;
			
			// run through the argument and remember all arguments we find (format: ${argumentName})
			for ( int j=0; j<length; j++ ) {
				if ( start < 0 ) {
					// we haven't found the start of an argument yet
					if ( value[j] == '$' && (j+1) < length && value[j+1] == '{' ) {
						// we found an argument
						start = j+2;
						j++;
					}
				} else if ( value[j] == '}' ) {
					// we found the end of the argument, so save it
					size_t argLength = j-start;
					char *argName = (char*)malloc(sizeof(char)*(argLength+1));
					sprintf(argName, "%.*s", (int)argLength, value + start);
					allMissingArguments.push_back(argName);
					start = -1;
				}
			}
		}
	}
	
	// only keep every argument once
	vector<char*> uniqueArguments;
	for ( vector<char*>::size_type t = 0; t < allMissingArguments.size(); t++ ) {
		char *missingArgument = allMissingArguments[t];
		
		bool found = false;
		for ( vector<char*>::size_type i = 0; i < uniqueArguments.size(); i++ ) {
			if ( ! strcmp(uniqueArguments[i], missingArgument) ) {
				found = true;
				break;
			}
		}
		
		if ( ! found ) {
			uniqueArguments.push_back(missingArgument);
		}
	}
	
	return uniqueArguments;
}

char* RSJobParser::leftTrimString(char *s)
{
    while(isspace(*s)) s++;
    return s;
}

char* RSJobParser::rightTrimString(char *s)
{
    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}

char* RSJobParser::trimString(char *s)
{
    return this->rightTrimString(this->leftTrimString(s));
}

// from: http://stackoverflow.com/questions/779875/what-is-the-function-to-replace-string-in-c
char* RSJobParser::replaceString(char *orig, char *rep, char *with)
{
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep
    int len_with; // length of with
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    if (!orig)
        return NULL;
    if (!rep)
        sprintf(rep, "");
    len_rep = strlen(rep);
    if (!with)
		sprintf(with, "");
    len_with = strlen(with);

    ins = orig;
    for (count = 0; tmp = strstr(ins, rep); ++count) {
        ins = tmp + len_rep;
    }

    // first time through the loop, all the variable are set correctly
    // from here on,
    //    tmp points to the end of the result string
    //    ins points to the next occurrence of rep in orig
    //    orig points to the remainder of orig after "end of rep"
    tmp = result = (char*)malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);
    return result;
}
