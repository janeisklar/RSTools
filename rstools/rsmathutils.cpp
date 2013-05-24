/*******************************************************************
 *
 * rsmathutils.cpp
 *
 * 
 * André Hoffmann
 *******************************************************************/

#include <RInside.h>
#include "rsmathutils.h"

int testMath()
{
    const char * const argv = "";
    RInside R(0, &argv, true, true, false);              // create an embedded R instance

    R["txt"] = "Hello, world!\n";	// assign a char* (string) to 'txt'

    R.parseEvalQ("cat(txt)");           // eval the init string, ignoring any returns

    return 3;
}
