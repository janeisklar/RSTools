#include "rsapplytransformation_common.h"
#include "rsapplytransformation_ui.h"

int main(int argc, char * argv[])
{
    // Parse run arguments
    rsApplyTransformationParameters *p = rsApplyTransformationParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsApplyTransformationInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsApplyTransformationRun(p);
    }
    
    BOOL execSuccessful = p->parametersValid;

    // Free memory
    rsApplyTransformationDestroy(p);

    return execSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
