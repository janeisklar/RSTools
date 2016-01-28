#include "rssmoothing_common.h"
#include "rssmoothing_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsSmoothingParameters *p = rsSmoothingParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsSmoothingInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsSmoothingRun(p);
    }
    
    BOOL execSuccessful = p->parametersValid;

    // Free memory
    rsSmoothingDestroy(p);

    return execSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
