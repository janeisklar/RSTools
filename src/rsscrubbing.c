#include "rsscrubbing_common.h"
#include "rsscrubbing_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsScrubbingParameters *p = rsScrubbingParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsScrubbingInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsScrubbingRun(p);
    }
    
    BOOL execSuccessful = p->parametersValid;

    // Free memory
    rsScrubbingDestroy(p);

    return execSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
