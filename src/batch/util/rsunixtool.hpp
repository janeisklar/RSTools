#ifndef rstools_rsbatch_util_unixtool_h
#define rstools_rsbatch_util_unixtool_h

#include "rstool.hpp"
#include "rsunixtask.hpp"
#include <vector>

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
        virtual bool _setupTempDir();
        virtual bool _prepareRun();
        virtual bool _prepareStream();
        virtual void _run();
        virtual void _moveOutputIfNecessary();
        virtual void _finalizeRun();
        virtual bool _createStream(const char* path);
        virtual RSUnixTask* getUnixTask();

        bool executionSuccessful;
        char* streamName;
        char* streamTarget;
        char* tmpDirPath;
        char* executionCmd;
        char* executionCmdPrint;

        nifti_image *inputNifti;

    private:
        bool _interceptOutput(FILE *stream);
};

}}} // namespace rstools::batch::util

#endif
