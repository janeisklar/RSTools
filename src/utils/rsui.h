#ifndef rstools_utils_ui_h
#define rstools_utils_ui_h

#include "src/rscommon.h"
#include "rsstring.h"
#include "glib.h"

#ifdef __cplusplus
extern "C" {
#endif

static const short RS_UI_GROUP_MAIN     = 1;
static const short RS_UI_GROUP_EXTENDED = 2;

typedef struct {
    
    // GOptionEntry varibles
    gchar          *name;
    gchar          shorthand;
    gint           flags;

    GOptionArg     type;
    gpointer       storage;
  
    gchar*         cli_description;
    gchar*         cli_arg_description;
    
    // Additional variables
    gchar*         gui_description;
    unsigned short group;
    
} rsUIOption;

typedef struct {
    gchar      *description;
    rsUIOption **options;
    size_t     nOptions;
} rsUIInterface;

rsUIInterface*  rsUINewInterface();
void            rsUIDestroyInterface(rsUIInterface* I);

rsUIOption*     rsUINewOption();
void            rsUIDestroyOption(rsUIOption* o);

void            rsUIAddOption(rsUIInterface* I, rsUIOption* o);

GOptionContext* rsUICreateCLI(rsUIInterface* I);
BOOL            rsUIParse(rsUIInterface* I, int argc, char * argv[]);

char*           rsUIPrepareHelpText(const char* text, const size_t lineWidth);

#ifdef __cplusplus
}
#endif

#endif
