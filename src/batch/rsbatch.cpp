#include <stdio.h>
#include "rsbatch_common.hpp"
#include "rsbatch_ui.hpp"

int main( int argc, char* argv[] )
{
	// Parse run arguments
	rsBatchParameters *p = rsBatchParseParams(argc, argv);
	
	if( argc < 2 ) {
		rsBatchPrintHelp(p);
		return 0;
	}
	
	// If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
		rsBatchInit(p);
    }
	
	// If everything went well start the main processing
	if ( p->parametersValid ) {
		rsBatchRun(p);
    }

	// Free memory
    rsBatchDestroy(p);

	return 0;
}
