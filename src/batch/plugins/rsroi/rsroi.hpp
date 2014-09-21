#ifndef rstools_rsbatch_plugins_rsroi_rsroi_h
#define rstools_rsbatch_plugins_rsroi_rsroi_h

#include "rscommon.h"
#include "batch/util/plugin.hpp"
#include "tool/roi.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsroi {

class RSRoi : Plugin {

    public:
        RSRoi();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createRoiTool();
        static RSTask* createRoiTask();
        
    protected:
        rsToolRegistration* createRoiToolRegistration();
        rsXSDExtension* createRoiToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rsroi

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
