#include "rsresample_common.h"
#include "rsresample_ui.h"

int main(int argc, char * argv[])
{    
    // Parse run arguments
    rsResampleParameters *p = rsResampleParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsResampleInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsResampleRun(p);
    }
    
    BOOL execSuccessful = p->parametersValid;

    // Free memory
    rsResampleDestroy(p);

    return execSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
