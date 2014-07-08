#ifndef rstools_rsjob_hpp
#define rstools_rsjob_hpp

#include <vector>
#include "rstask.hpp"

using namespace std;

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

#endif
