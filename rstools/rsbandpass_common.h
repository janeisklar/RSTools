//
//  rsbandpass_common.h
//  rstools
//
//  Created by André Hoffmann on 6/28/13.
//
//

#ifndef rstools_rsbandpass_common_h
#define rstools_rsbandpass_common_h

#include <stdio.h>
#include <nifti1.h>
#include <fslio.h>
#include "rsniftiutils.h"
#include "rsmathutils.h"

struct rsBandpassParameters {
    char *inputpath;
	char *maskpath;
    char *savemaskpath;
    char *saveFilteredPath;
    char *saveAttenuationPath;
	
	short xDim;
    short yDim;
    short zDim;
    short vDim;
    
	short pixtype;
	size_t dt;
    
    float inter;
    float slope;
    
    double f1;
    double f2;
    double sampling_rate;
    
    BOOL verbose;
    BOOL parametersValid;
    
    FSLIO *fslio;
    double ***mask;

    int threads;
    size_t wordsize;
    
    struct rsFFTFilterParams fftParams;
    double rolloff;
    int rolloff_method;
    long paddedT;
};

struct rsBandpassParameters rsBandpassLoadParams(int argc, char * argv[]);
struct rsBandpassParameters rsBandpassInitParameters();
void rsBandpassPrintHelp();
void testFFTFilter();

#endif
