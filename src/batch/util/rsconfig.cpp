#include "rsconfig.hpp"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

RSConfig::RSConfig()
{
    _loadConfig();
}

RSConfig::~RSConfig()
{}

void RSConfig::_loadConfig()
{
    const char *confPath = CONFIG_PATH"/rstools.conf";
    ifstream c(confPath);
    
    int n = 1;
    string line;
    while (getline(c, line)) {
        
        if ( line[0] == '#' || line.find("=")==string::npos ) {
            continue;
        }
        
        size_t delimiter = line.find("=");
        
        if ( delimiter == string::npos ) {
            fprintf(stderr, "Error while parsing config '%s' on line %d:\n%s\n", confPath, n, line.c_str());
            exit(EXIT_FAILURE);
        }
        
        string key = line.substr(0, delimiter);
        string value = line.substr(delimiter+1);
        
        rsArgument *arg = (rsArgument*)rsMalloc(sizeof(rsArgument));
        arg->key = rsString(key.c_str());
        arg->value = rsString(value.c_str());
        arguments.push_back(arg);
        
        n++;
    }
}

void RSConfig::reload()
{
    arguments.clear();
    _loadConfig();
}

void RSConfig::setArguments(vector<rsArgument*> newArguments)
{
    arguments.clear();

    for ( vector<rsArgument*>::size_type i = 0; i != newArguments.size(); i++ ) {
        arguments.push_back(newArguments[i]);
    }
}

vector<rsArgument*> RSConfig::getArguments()
{
    return this->arguments;
}

rsArgument* RSConfig::getArgument(const char* name)
{
    for ( vector<rsArgument*>::size_type i = 0; i != arguments.size(); i++ ) {
        rsArgument *arg = arguments[i];
        
        if ( ! strcmp(arg->key, name) ) {
            return  arg;
        }
    }
    
    return NULL;
}

}}} // namespace rstools::batch::util
