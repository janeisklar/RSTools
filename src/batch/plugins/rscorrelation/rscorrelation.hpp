#ifndef rstools_rsbatch_plugins_rsbandpass_rscorrelation_h
#define rstools_rsbatch_plugins_rsbandpass_rscorrelation_h

#include "src/rscommon.h"
#include "src/batch/util/plugin.hpp"
#include "tool/correlation.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rscorrelation {

class RSCorrelation : Plugin {

    public:
        RSCorrelation();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createCorrelationTool();
        static RSTask* createCorrelationTask();
        
    protected:
        rsToolRegistration* createCorrelationToolRegistration();
        rsXSDExtension* createCorrelationToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rscorrelation

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
