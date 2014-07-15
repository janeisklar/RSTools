#ifndef rstools_rsbatch_plugins_rsbandpass_rsregression_h
#define rstools_rsbatch_plugins_rsbandpass_rsregression_h

#include "src/rscommon.h"
#include "src/batch/util/plugin.hpp"

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
        const char* getVersion();
};

}}}} // namespace rstools::batch::plugins::rsregression

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
