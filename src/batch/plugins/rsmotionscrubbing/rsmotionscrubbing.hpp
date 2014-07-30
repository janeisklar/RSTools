#ifndef rstools_rsbatch_plugins_rsmotionscrubbing_rsmotionscrubbing_h
#define rstools_rsbatch_plugins_rsmotionscrubbing_rsmotionscrubbing_h

#include "src/rscommon.h"
#include "src/batch/util/plugin.hpp"
#include "tool/motionscrubbing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsmotionscrubbing {

class RSMotionScrubbing : Plugin {

    public:
        RSMotionScrubbing();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createMotionScrubbingTool();
        static RSTask* createMotionScrubbingTask();
        
    protected:
        rsToolRegistration* createMotionScrubbingToolRegistration();
        rsXSDExtension* createMotionScrubbingToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rsmotionscrubbing

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
