#ifndef rstools_rsbatch_util_rsunixtask_hpp
#define rstools_rsbatch_util_rsunixtask_hpp

#include <nifti/headerinfo.h>
#include "rstask.hpp"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

class RSUnixTask: public RSTask {
    
    public:
        
        RSUnixTask(const char*, const char*);

        virtual char* getCmd(bool asExecuted) = 0;

        virtual bool hasInputNiftiHeaderInformation();
        virtual rsNiftiExtendedHeaderInformation* getInputNiftiHeaderInformation();
        virtual void setInputNiftiHeaderInformation(rsNiftiExtendedHeaderInformation* info);

        virtual char *getTempDirectoryPath();
        virtual void setTempDirectoryPath(char *path);

    protected:
        rsNiftiExtendedHeaderInformation* inputNiftiHeaderInformation;
        char *tempDirectoryPath;
};

}}} // namespace rstools::batch::util

#endif
