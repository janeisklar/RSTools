#include "rsui.h"

rsUIInterface* rsUINewInterface()
{
    rsUIInterface *I = rsMalloc(sizeof(rsUIInterface));
    
    I->nOptions    = 0;
    I->options     = (rsUIOption**)rsMalloc(sizeof(rsUIOption*));
    I->description = "Description: TBD";
    
    return I;
}

void rsUIDestroyInterface(rsUIInterface* I)
{
    for ( size_t i=0; i<I->nOptions; i++ ) {
        rsUIDestroyOption(I->options[i]);
    }
    rsFree(I);
}

rsUIOption* rsUINewOption()
{
    rsUIOption *o = (rsUIOption*)rsMalloc(sizeof(rsUIOption));
    
    o->name                = NULL;
    o->shorthand           = 0;
    o->flags               = 0;
    o->type                = G_OPTION_ARG_NONE;
    o->storage             = NULL;
    o->cli_description     = NULL;
    o->cli_arg_description = NULL;
    o->gui_description     = NULL;
    o->group               = RS_UI_GROUP_MAIN;
    o->showInGUI           = TRUE;
    
    return o;
}

void rsUIDestroyOption(rsUIOption* o)
{
    rsFree(o->gui_description);
    rsFree(o);
}

void rsUIAddOption(rsUIInterface* I, rsUIOption* o)
{
    I->nOptions++;
    I->options = (rsUIOption**)realloc(I->options, I->nOptions * sizeof(rsUIOption*));
    I->options[I->nOptions-1] = o;
}

GOptionContext* rsUICreateCLI(rsUIInterface* I)
{
    // create a new context using the description
    char *lineFeed = "\n";
    char *description = rsUIPrepareHelpText(I->description, 60, lineFeed);
    char *oldDescription = description;
    
    description = rsStringConcat("", description, NULL);
    rsFree(oldDescription);
    
    GOptionContext *context = g_option_context_new(description);
    rsFree(description);
    
    // add the current version number to the header output of the help
    g_option_context_set_summary(context, RSTOOLS_VERSION_LABEL);
    
    // count the number of options in the main and extended group
    size_t nMainOptions = 0, nExtendedOptions = 0;
    
    for ( size_t i=0; i<I->nOptions; i++ ) {
        if (I->options[i]->group == RS_UI_GROUP_EXTENDED ) {
            nExtendedOptions++;
        } else {
            nMainOptions++;
        }
    }
    
    // create option structs for all options
    GOptionEntry *mainEntries     = rsMalloc(sizeof(GOptionEntry)*(nMainOptions+1));
    GOptionEntry *extendedEntries = rsMalloc(sizeof(GOptionEntry)*(nExtendedOptions+1));
    
    GOptionEntry *iMainEntries     = mainEntries;     // iterators
    GOptionEntry *iExtendedEntries = extendedEntries;
    
    for ( size_t i=0; i<I->nOptions; i++ ) {
        rsUIOption* o = I->options[i];

        // limit line width so that the description is not too wide
        char *lineFeed = "\n                                   ";
        description = rsUIPrepareHelpText(o->cli_description, 45, lineFeed);
        
        // append a new line for better readability
        oldDescription = description;
        description = rsStringConcat(description, "\n", NULL);
        rsFree(oldDescription);

        GOptionEntry entry = {
            o->name,
            o->shorthand,
            o->flags,
            o->type,
            o->storage,
            description, 
            o->cli_arg_description
        };
        
        GOptionEntry *entries = iMainEntries;
        
        if (o->group == RS_UI_GROUP_EXTENDED ) {
            entries = iExtendedEntries;
        }
        
        memcpy(entries, &entry, sizeof(GOptionEntry));
            
        if (o->group == RS_UI_GROUP_MAIN ) {
            iMainEntries++;
        } else {
            iExtendedEntries++;
        }
    }
    
    // terminate the array of structs with a NULL element as required by glib
    GOptionEntry entry = {NULL};
    memcpy(iMainEntries,     &entry, sizeof(GOptionEntry));
    memcpy(iExtendedEntries, &entry, sizeof(GOptionEntry));

    // add main options
    g_option_context_add_main_entries(context, mainEntries, GETTEXT_PACKAGE);
    
    if ( nExtendedOptions < 1 ) {
        return context;
    }
    
    // add extended options if there are any
    GOptionGroup *g = g_option_group_new("extended", "Extended options", "Additional options that are rarely going to be used", NULL, NULL);
    g_option_context_add_group(context, g);
    g_option_group_add_entries(g, extendedEntries);
    
    return context;
}

BOOL rsUIParse(rsUIInterface* I, int argc, char * argv[])
{
    GOptionContext* context = rsUICreateCLI(I);
    GError*         error   = NULL;
    
    if ( argc < 2 ) {
        fprintf(stdout, "%s\n", g_option_context_get_help(context, TRUE, NULL));
        return FALSE;
    }
    
    if ( ! g_option_context_parse(context, &argc, &argv, &error) ) {
        fprintf(stderr, "option parsing failed: %s\n", error->message);
        return FALSE;
    }
    
    return TRUE;
}

char* rsUIPrepareHelpText(const char* text, const size_t lineWidth, char *lineGlue)
{
    char **textArray;
    size_t textArraySize;
    char *helpText = (char*)rsMalloc(sizeof(char)); helpText[0] = '\0';
    rsStringWordWrap(text, &textArray, &textArraySize, lineWidth);
    
    for ( size_t j=0; j< textArraySize; j++) {
        char *oldText = helpText;
        char *glue = (j==0) ? "" : lineGlue;
        helpText = rsStringConcat(helpText, glue, textArray[j], NULL);
        rsFree(oldText);
    }
    
    return helpText;
}
