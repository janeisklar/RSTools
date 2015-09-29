#include "rsinfo_common.h"
#include "rsinfo_ui.h"

int main(int argc, char * argv[]) {

    // Parse run arguments
    rsInfoParameters *p = rsInfoParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsInfoInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsInfoRun(p);
    }

    // Free memory
    rsInfoDestroy(p);

    return 0;
}
