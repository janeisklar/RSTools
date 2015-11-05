#include "rsdeoblique_common.h"
#include "rsdeoblique_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsDeobliqueParameters *p = rsDeobliqueParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsDeobliqueInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsDeobliqueRun(p);
    }

    // Free memory
    rsDeobliqueDestroy(p);

    return 0;
}
