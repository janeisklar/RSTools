#ifndef rstools_rsbatch_util_unixtool_h
#define rstools_rsbatch_util_unixtool_h

#include "rstool.hpp"
#include "rsunixtask.hpp"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

class RSUnixTool : public RSTool {

    public:    
        virtual void printCallString(FILE *stream);
        virtual bool isEverythingFine();
        
    protected:
        virtual void _parseParams(int argc, char * argv[]);
        virtual void _run();
        virtual RSUnixTask* getUnixTask();
        
        bool executionSuccessful;
};

}}} // namespace rstools::batch::util

#endif
