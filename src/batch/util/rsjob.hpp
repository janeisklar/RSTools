#ifndef rstools_rsbatch_util_rsjob_hpp
#define rstools_rsbatch_util_rsjob_hpp

#include <vector>
#include "rstask.hpp"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

class RSJob {
    protected:
        char *jobfile;
        char *description;
        vector<RSTask*> tasks;
        vector<rsArgument*> arguments;
        
    public:
        RSJob(char *jobfile);
        ~RSJob();
        
        void setJobfile(char*);
        char* getJobfile();
        void setDescription(char*);
        char* getDescription();
        void addTask(RSTask*);
        vector<RSTask*> getTasks();
        void addArgument(rsArgument*);
        vector<rsArgument*> getArguments();
};

}}} // namespace rstools::batch::util

#endif
