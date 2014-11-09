#include "rsjobparser.hpp"
#include "rsconfig.hpp"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

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
    
    const char *schemaFile = mergePluginXSDExtensions();
    const char schemaLocationTemplate[] = "http://www.fmri.at/rstools %s";
    char schemaLocation[strlen(schemaLocationTemplate)+strlen(schemaFile)+1];
    sprintf(&schemaLocation[0], schemaLocationTemplate, schemaFile);
    
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
    
    unlink(schemaFile);

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
        arg->key   = rsTrimString(XMLString::transcode(current_node->getAttributes()->getNamedItem(XMLString::transcode("name"))->getNodeValue()));
        arg->value = rsTrimString(XMLString::transcode(current_node->getFirstChild()->getNodeValue()));
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
        char *taskName = XMLString::transcode(current_node->getNodeName());
        RSTask *task = RSTask::taskFactory(taskName);
        task->parseTaskFromXml(walker, current_node);
        this->job->addTask(task);
    } while ( current_node != NULL && ! strcmp("tasks", XMLString::transcode(current_node->getParentNode()->getNodeName())) );
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
    
    vector<rsArgument*> configArguments = RSConfig::getInstance().getArguments();
    
    // merge job arguments and config arguments
    for ( int u=0; u<configArguments.size(); u++ ) {
        rsArgument* configArgument = configArguments[u];
        bool found = false;

        // check if config argument already exists in the list of job arguments
        for ( vector<rsArgument*>::size_type i = 0; i < jobArguments.size(); i++ ) {
            rsArgument *arg = jobArguments[i];

            if ( strcmp(arg->key, configArgument->key) ) {
                continue;
            }

            // if it already exists use the existing value
            found = true;
            break;
        }

        // if it doesn't exist yet add it to the job arguments
        if ( ! found ) {
            rsArgument* newArgument = (rsArgument*)malloc(sizeof(rsArgument));
            newArgument->key   = (char*)malloc((strlen(configArgument->key)  +(size_t)1)*sizeof(char));
            newArgument->value = (char*)malloc((strlen(configArgument->value)+(size_t)1)*sizeof(char));
            sprintf(newArgument->key,   "%s", configArgument->key);
            sprintf(newArgument->value, "%s", configArgument->value);
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
    
    RSConfig::getInstance().setArguments(jobArguments);

    // replace all job argument placeholders within the task arguments
    vector<RSTask*> tasks = this->job->getTasks();
    for ( vector<RSTask*>::size_type t = 0; t < tasks.size(); t++ ) {
        RSTask *task = tasks[t];
        task->fillInJobArguments(job, this);
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

// from: http://stackoverflow.com/questions/779875/what-is-the-function-to-replace-string-in-c
char* RSJobParser::replaceString(const char *orig, const char *rep, const char *with)
{
    char *result;  // the return string
    char *ins;     // the next insert point
    char *tmp;     // varies
    char *input;   // copy of orig
    int len_rep;   // length of rep
    int len_with;  // length of with
    int len_front; // distance between rep and end of last rep
    int count;     // number of replacements

    if (!orig)
        return NULL;
    
    // copy the original string so it is not affected
    input = (char*)malloc(sizeof(char)*(strlen(orig)+1));
    sprintf(input, "%s", orig);
    
    len_rep = strlen(rep);
    len_with = strlen(with);

    ins = input;
    
    for (count = 0; (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }

    // first time through the loop, all the variable are set correctly
    // from here on,
    //    tmp points to the end of the result string
    //    ins points to the next occurrence of rep in input
    //    input points to the remainder of orig after "end of rep"
    tmp = result = (char*)malloc(strlen(input) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    while (count--) {
        ins = strstr(input, rep);
        len_front = ins - input;
        tmp = strncpy(tmp, input, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        input += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, input);
    return result;
}

/*
 * Merges all known extension of the XSD schema
 * into a single temporary file and returns its
 * file path
 */
char* RSJobParser::mergePluginXSDExtensions()
{
    vector<rsXSDExtension*> extensions = RSTool::getXSDExtensions();
    stringstream extensionStream;
    stringstream typeStream;
    
    // read in and merge all extensions
    for ( vector<rsXSDExtension*>::size_type t = 0; t < extensions.size(); t++ ) {
        rsXSDExtension *extension = extensions[t];
        const char* filePath = extension->file;
        
        ifstream extfile;
        extfile.open(filePath, ifstream::in | ios::binary);
        
        if ( extfile.is_open() ) {
            extensionStream << extfile.rdbuf() << endl;            
            extfile.close();
        } else {
            throw runtime_error("Could not read a plugin's job validation extension file!");
        }
    }
    
    string extensionReplacement = extensionStream.str();
    
    // create list of registered types
    for ( vector<rsXSDExtension*>::size_type t = 0; t < extensions.size(); t++ ) {
        rsXSDExtension *extension = extensions[t];
        typeStream << "                <xs:element name=\"" << extension->name << "\" type=\"" << extension->type << "\"/>" << endl;
    }
    
    string typeReplacement = typeStream.str();
    
    // read in main xsd-schema definiton
    const char* xsdFilePath = RSTOOLS_DATA_DIR "/" PACKAGE "/jobs/job.xsd";
    stringstream schemaStream;
    ifstream xsdfile;
    xsdfile.open(xsdFilePath, ifstream::in | ios::binary);

    if ( xsdfile.is_open() ) {
        schemaStream << xsdfile.rdbuf() << endl;
        xsdfile.close();
    } else {
        throw runtime_error("Could not read job valdation file!");
    }
    
    string xsdSchema = schemaStream.str();
        
    // replace both the typelist and the extensions in the xsdSchema
    replaceAll(xsdSchema, string("<!-- %%TYPELIST%% -->"  ),      typeReplacement);
    replaceAll(xsdSchema, string("<!-- %%EXTENSIONS%% -->"), extensionReplacement);
    
    // write result to a temporary file
    char* tmpFileName = (char*)malloc(sizeof(char)*255);
    sprintf(tmpFileName, "%s", "/tmp/job.xsd-XXXXXXX");
    int filedes = mkstemp(tmpFileName);
    
    if ( filedes < 1 ) {
        throw runtime_error("Could not create a temporary file to merge the job-validation files for all plugins.");
    }
    
    if ( -1 == write(filedes, xsdSchema.c_str(), strlen(xsdSchema.c_str())) ) {
        throw runtime_error("Could not write to temporary file");
    }
    
    // return path to the merged xsd-file
    return tmpFileName;
}

}}} // namespace rstools::batch::util
