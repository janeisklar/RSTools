#include <src/utils/rsui.h>
#include "rsmaskborderdistance_common.h"
#include "rsmaskborderdistance_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsMaskBorderDistanceParameters *p = rsMaskBorderDistanceParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsMaskBorderDistanceInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsMaskBorderDistanceRun(p);
    }
    
    BOOL execSuccessful = p->parametersValid;

    // Free memory
    rsMaskBorderDistanceDestroy(p);

    return execSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
