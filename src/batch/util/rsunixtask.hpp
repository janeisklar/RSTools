#ifndef rstools_rsbatch_util_rsunixtask_hpp
#define rstools_rsbatch_util_rsunixtask_hpp

#include "rstask.hpp"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

class RSUnixTask: public RSTask {
    
    public:
        
        RSUnixTask(const char*, const char*);
        
        virtual char* getCmd() = 0;
    
    protected:
        
};

}}} // namespace rstools::batch::util

#endif
