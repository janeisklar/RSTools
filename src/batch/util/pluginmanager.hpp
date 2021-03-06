#ifndef rstools_rsbatch_util_pluginmanager_h
#define rstools_rsbatch_util_pluginmanager_h

#include <cstdio>
#include <glib.h>
#include <gmodule.h>
#include <glob.h>
#include <vector>
#include <string>
#include <stdexcept>
#include "plugin.hpp"
#include "rscommon.h"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

class PluginManager {
        
    public:
        static PluginManager& getInstance();
        void loadPlugins();
        vector<string> getPluginLibraries();
        void loadLibrary(const char* library);
        void registerPlugin(Plugin* plugin);
        
    private:
        PluginManager()
	{
		pluginsLoaded = false;
	}; 
        PluginManager(PluginManager const&);
        void operator=(PluginManager const&);

    protected:
        bool pluginsLoaded;
        vector<Plugin*> plugins;
        
        
};

}}} // namespace rstools::batch::util

#endif
