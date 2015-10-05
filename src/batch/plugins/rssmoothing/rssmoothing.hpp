#ifndef rstools_rsbatch_plugins_rssmoothing_rssmoothing_h
#define rstools_rsbatch_plugins_rssmoothing_rssmoothing_h

#include "rscommon.h"
#include "batch/util/plugin.hpp"
#include "tool/smoothing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rssmoothing {

class RSSmoothing : Plugin {

    public:
        RSSmoothing();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createSmoothingTool();
        static RSTask* createSmoothingTask();
        
    protected:
        rsToolRegistration* createSmoothingToolRegistration();
        rsXSDExtension* createSmoothingToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rssmoothing

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
