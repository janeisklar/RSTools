#include "rszeropadding_common.h"
#include "rszeropadding_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsZeropaddingParameters *p = rsZeropaddingParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsZeropaddingInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsZeropaddingRun(p);
    }
    
    BOOL execSuccessful = p->parametersValid;

    // Free memory
    rsZeropaddingDestroy(p);

    return execSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
