#include "rscorrelation_common.h"
#include "rscorrelation_ui.h"

int main(int argc, char * argv[]) {

    // Parse run arguments
    rsCorrelationParameters *p = rsCorrelationParseParams(argc, argv);
    
    // If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
        rsCorrelationInit(p);
    }
    
    // If everything went well start the main processing
    if ( p->parametersValid ) {
        rsCorrelationRun(p);
    }

    // Free memory
    rsCorrelationDestroy(p);

    return 0;
}
