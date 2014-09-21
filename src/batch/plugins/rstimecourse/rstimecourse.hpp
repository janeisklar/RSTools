#ifndef rstools_rsbatch_plugins_rstimecourse_rstimecourse_h
#define rstools_rsbatch_plugins_rstimecourse_rstimecourse_h

#include "rscommon.h"
#include "batch/util/plugin.hpp"
#include "tool/timecourse.hpp"

using namespace rstools::batch::util;

namespace rstools {
namespace batch {
namespace plugins {
namespace rstimecourse {

class RSTimecourse : Plugin {

    public:
        RSTimecourse();
        void  registerPlugin();
        const char* getName();
        const char* getCode();
        const char* getVersion();

        static RSTool* createTimecourseTool();
        static RSTask* createTimecourseTask();
        
    protected:
        rsToolRegistration* createTimecourseToolRegistration();
        rsXSDExtension* createTimecourseToolXSDExtension();
};

}}}} // namespace rstools::batch::plugins::rstimecourse

extern "C" {
    G_MODULE_EXPORT Plugin* rsGetPlugin(void);
}

#endif
