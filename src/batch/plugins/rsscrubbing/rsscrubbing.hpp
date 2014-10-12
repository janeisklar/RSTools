#ifndef rstools_rsbatch_plugins_rsscrubbing_rsscrubbing_h
#define rstools_rsbatch_plugins_rsscrubbing_rsscrubbing_h

#include "rscommon.h"
#include "batch/util/plugin.hpp"
#include "tool/scrubbing.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsscrubbing {

class RSScrubbing : Plugin {

    public:
        RSScrubbing();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createScrubbingTool();
        static RSTask* createScrubbingTask();
        
    protected:
        rsToolRegistration* createScrubbingToolRegistration();
        rsXSDExtension* createScrubbingToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rsscrubbing

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
