#ifndef rstools_rsbatch_plugins_rsapplytransformation_rsapplytransformation_h
#define rstools_rsbatch_plugins_rsapplytransformation_rsapplytransformation_h

#include "rscommon.h"
#include "batch/util/plugin.hpp"
#include "tool/applytransformation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsapplytransformation {

class RSApplyTransformation : Plugin {

    public:
        RSApplyTransformation();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createApplyTransformationTool();
        static RSTask* createApplyTransformationTask();
        
    protected:
        rsToolRegistration* createApplyTransformationToolRegistration();
        rsXSDExtension* createApplyTransformationToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rsapplytransformation

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
