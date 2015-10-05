#ifndef rstools_rsbatch_plugins_rsorientation_rsorientation_h
#define rstools_rsbatch_plugins_rsorientation_rsorientation_h

#include "rscommon.h"
#include "batch/util/plugin.hpp"
#include "tool/orientation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsorientation {

class RSOrientation : Plugin {

    public:
        RSOrientation();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createOrientationTool();
        static RSTask* createOrientationTask();
        
    protected:
        rsToolRegistration* createOrientationToolRegistration();
        rsXSDExtension* createOrientationToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rsorientation

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
