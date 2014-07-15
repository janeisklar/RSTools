#include "pluginmanager.hpp"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

void PluginManager::loadPlugins()
{
    if ( pluginsLoaded ) {
        return;
    }
    
    vector<string> libs = getPluginLibraries();
    
    for(vector<string>::iterator it = libs.begin(); it != libs.end(); ++it) {
        const char* lib = (*it).c_str();
        loadLibrary(lib);
    }
    
    pluginsLoaded = true;
}

void PluginManager::loadLibrary(const char* library)
{
    fprintf(stdout, "Loading Plugin: %s\n", library);
    
    GModule *module = NULL;
    rsGetPluginFunc getPluginFunc = NULL;

    /* Check whether glib is compiled with plugin support */
    if(g_module_supported() == FALSE) {
        throw runtime_error("GLib is not compiled with module support which is necessary to load the tools");
    }

    module = g_module_open(library, G_MODULE_BIND_MASK);

    /* Check whether the plugin was loaded successfully */
    if(module == NULL) {
        throw runtime_error( (string("Library '") + string(library) + string("' could not be opened, because of the following error: ") + string(g_module_error())).c_str() );
    }

    /* Load the symbol and assign it to our function pointer */
    if(g_module_symbol(module, "rsGetPlugin", (gpointer *) &getPluginFunc) == FALSE) {
        throw runtime_error( (string("Library '") + string(library) + string("' could not be loaded, due to the following error: ") + string(g_module_error()) ).c_str() );
    }

    /* acquire the plugin */
    Plugin *plugin = getPluginFunc();
    registerPlugin(plugin);

    //if(g_module_close(module) == FALSE) {
    //    throw runtime_error( (string("Library '") + string(library) + string("' could not be closed, due to the following error: ") + string(g_module_error()) ).c_str() );
    //}
}

vector<string> PluginManager::getPluginLibraries()
{
    vector<string> libraries;
    
    const char* pattern    = "/lib*Plugin.so";
    char* searchPath = (char*)malloc( sizeof(char) * (strlen(RSTOOLS_PLUGINS_DIR) + strlen(pattern) + 1) );
    sprintf(searchPath, "%s%s", RSTOOLS_PLUGINS_DIR, pattern);
        
    glob_t glob_result;
    glob(searchPath, GLOB_TILDE, NULL, &glob_result);
    
    for (unsigned int i=0; i<glob_result.gl_pathc; ++i) {
        libraries.push_back(string(glob_result.gl_pathv[i]));
    }
    
    globfree(&glob_result);
    free(searchPath);
    
    return libraries;
}

void PluginManager::registerPlugin(Plugin* plugin)
{
    // ensure that we haven't already loaded the plugin
    for(vector<Plugin*>::iterator it = plugins.begin(); it != plugins.end(); ++it) {
        Plugin *p = *it;
        if ( ! strcmp(p->getName(), plugin->getName()) ) {
            throw runtime_error( (string("Plugin with name '") + string(plugin->getName()) + string("' was loaded multiple times!")) );
        }
    }
    
    plugins.push_back(plugin);
    plugin->registerPlugin();
}

PluginManager& PluginManager::getInstance()
{
    static PluginManager instance;
    return instance;
}

}}} // namespace rstools::batch::util
