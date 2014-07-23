#ifndef rstools_rsbatch_plugins_unix_unix_h
#define rstools_rsbatch_plugins_unix_unix_h

#include "src/rscommon.h"
#include "src/batch/util/plugin.hpp"
#include "tool/unix.hpp"
#include "task/unix.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace unix {

class Unix : Plugin {

    public:
        Unix();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createUnixTool();
        static RSTask* createUnixTask();
        
    protected:
        rsToolRegistration* createUnixToolRegistration();
        rsXSDExtension* createUnixToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::unix

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
