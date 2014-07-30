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
    char *name;
    char *description;
} rsUIOptionValue;

typedef struct {
    
    // GOptionEntry varibles
    gchar*            name;
    gchar             shorthand;
    gint              flags;

    GOptionArg        type;
    gpointer          storage;
  
    gchar*            cli_description;
    gchar*            cli_arg_description;
    
    // Additional variables
    gchar*            gui_description; // used in the GUI instead of cli_description if present
    unsigned short    group;           // whether this option is listed in the main or extended group (see RS_UI_GROUP_*)
    BOOL              showInGUI;       // whether the option should be visible in the GUI
    rsUIOptionValue** allowedValues;   // (NULL-terminated) list of values that this option accepts
    char*             defaultValue;
    
} rsUIOption;

typedef struct {
    char*        description;
    char*        gui_description;
    rsUIOption** options;
    size_t       nOptions;
    size_t       helpIndent;
} rsUIInterface;

rsUIInterface*  rsUINewInterface();
void            rsUIDestroyInterface(rsUIInterface* I);

rsUIOption*     rsUINewOption();
void            rsUIDestroyOption(rsUIOption* o);

void            rsUIAddOption(rsUIInterface* I, rsUIOption* o);
void            rsUISetOptionValues(rsUIOption* o, rsUIOptionValue values[]);

GOptionContext* rsUICreateCLI(rsUIInterface* I, void *userData);
BOOL            rsUIParse(rsUIInterface* I, int argc, char * argv[], void *userData);

char*           rsUIPrepareHelpText(const char* text, const size_t lineWidth, char *lineGlue);

#ifdef __cplusplus
}
#endif

#endif
