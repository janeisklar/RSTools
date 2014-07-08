#include "rsmotionscrubbing_common.h"
#include "rsmotionscrubbing_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsMotionScrubbingParameters *p = rsMotionScrubbingParseParams(argc, argv);
    
    if( argc < 2 ) {
        rsMotionScrubbingPrintHelp(p);
        return 0;
    }
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsMotionScrubbingInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsMotionScrubbingRun(p);
    }

    // Free memory
    rsMotionScrubbingDestroy(p);

    return 0;
}
