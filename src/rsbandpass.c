#include "rsbandpass_common.h"
#include "rsbandpass_ui.h"

int main(int argc, char * argv[])
{
	// Parse run arguments
	rsBandpassParameters *p = rsBandpassParseParams(argc, argv);
	
	if( argc < 2 ) {
		rsBandpassPrintHelp(p);
		return 0;
	}
	
	// If arguments are valid, initialize niftis, etc.
    if ( p->parametersValid ) {
		rsBandpassInit(p);
    }
	
	// If everything went well start the main processing
	if ( p->parametersValid ) {
		rsBandpassRun(p);
    }

	// Free memory
    rsBandpassDestroy(p);

	return 0;
}
