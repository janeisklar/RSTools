#ifndef rstools_rsbatch_util_plugin_h
#define rstools_rsbatch_util_plugin_h

#include <glib.h>
#include <gmodule.h>
#include "rstool.hpp"

namespace rstools {
namespace batch {
namespace util {

class Plugin {

    public:
        virtual const char* getName()        = 0;
        virtual const char* getVersion()     = 0;
        virtual void        registerPlugin() = 0;
        
};

typedef Plugin* (*rsGetPluginFunc) (void);

}}} // namespace rstools::batch::util

#endif
