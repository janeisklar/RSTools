#ifndef rstools_rsbatch_plugins_rsbandpass_rsregression_h
#define rstools_rsbatch_plugins_rsbandpass_rsregression_h

#include "rscommon.h"
#include "batch/util/plugin.hpp"
#include "tool/regression.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rsregression {

class RSRegression : Plugin {

    public:
        RSRegression();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createRegressionTool();
        static RSTask* createRegressionTask();
        
    protected:
        rsToolRegistration* createRegressionToolRegistration();
        rsXSDExtension* createRegressionToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rsregression

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
