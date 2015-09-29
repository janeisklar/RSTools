#include "rsorientation_common.h"
#include "rsorientation_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsOrientationParameters *p = rsOrientationParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsOrientationInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsOrientationRun(p);
    }

    // Free memory
    rsOrientationDestroy(p);

    return 0;
}
