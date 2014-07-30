#include "rsmotionscrubbing_common.h"
#include "rsmotionscrubbing_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsMotionScrubbingParameters *p = rsMotionScrubbingParseParams(argc, argv);
    
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
