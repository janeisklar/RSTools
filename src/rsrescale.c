#include "rsrescale_common.h"
#include "rsrescale_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsRescaleParameters *p = rsRescaleParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsRescaleInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsRescaleRun(p);
    }
    
    BOOL execSuccessful = p->parametersValid;

    // Free memory
    rsRescaleDestroy(p);

    return execSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
