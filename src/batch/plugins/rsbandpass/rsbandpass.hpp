#ifndef rstools_rsbatch_plugins_rsbandpass_rsbandpass_h
#define rstools_rsbatch_plugins_rsbandpass_rsbandpass_h

#include "src/rscommon.h"
#include "src/batch/util/plugin.hpp"
#include "tool/bandpass.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsbandpass {

class RSBandpass : Plugin {

    public:
        RSBandpass();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createBandpassTool();
        static RSTask* createBandpassTask();
        
    protected:
        rsToolRegistration* createToolRegistration();
};

}}}} // namespace rstools::batch::plugins::rsbandpass

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
