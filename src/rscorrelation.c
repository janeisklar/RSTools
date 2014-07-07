#include "rscorrelation_common.h"
#include "rscorrelation_ui.h"

int main(int argc, char * argv[]) {

	// Parse run arguments
    rsCorrelationParameters *p = rsCorrelationParseParams(argc, argv);
	
	if( argc < 2 ) {
		rsCorrelationPrintHelp(p);
		return 0;
	}
	
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
