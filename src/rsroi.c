#include "rsroi_common.h"
#include "rsroi_ui.h"

int main(int argc, char * argv[]) {

    // Parse run arguments
    rsRoiParameters *p = rsRoiParseParams(argc, argv);

    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsRoiInit(p);
    }

    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsRoiRun(p);
    }
    
    BOOL execSuccessful = p->parametersValid;

    // Free memory
    rsRoiDestroy(p);

    return execSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
