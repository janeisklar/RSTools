#include "rsregression_common.h"
#include "rsregression_ui.h"

int main(int argc, char * argv[]) {

    // Parse run arguments
    rsRegressionParameters *p = rsRegressionParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsRegressionInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsRegressionRun(p);
    }

    // Free memory
    rsRegressionDestroy(p);

    return p->parametersValid ? EXIT_SUCCESS : EXIT_FAILURE;
}
