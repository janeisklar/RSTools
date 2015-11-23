#include "rsbandpass_common.h"
#include "rsbandpass_ui.h"

int main(int argc, char * argv[])
{    
    // Parse run arguments
    rsBandpassParameters *p = rsBandpassParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsBandpassInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsBandpassRun(p);
    }

    // Free memory
    rsBandpassDestroy(p);

    return p->parametersValid ? EXIT_SUCCESS : EXIT_FAILURE;
}
