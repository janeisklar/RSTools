#ifndef rstools_rsbatch_util_rsconfig_hpp
#define rstools_rsbatch_util_rsconfig_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include "rstask.hpp"
#include "rscommon.h"
#include "utils/rsstring.h"

using namespace std;

namespace rstools {
namespace batch {
namespace util {
    
class RSConfig {
    
public:
    vector<rsArgument*> getArguments();
    static RSConfig& getInstance()
    {
        static RSConfig instance;
        return instance;
    }
    
protected:
    vector<rsArgument*> arguments;
    void _loadConfig();
        
private:
    RSConfig();
    ~RSConfig();
    RSConfig(RSConfig const&);
    void operator=(RSConfig const&);
};

}}} // namespace rstools::batch::util

#endif
