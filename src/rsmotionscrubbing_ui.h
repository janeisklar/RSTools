#ifndef rstools_motionscrubbing_ui_h
#define rstools_motionscrubbing_ui_h

#include <stdio.h>
#include <nifti1.h>
#include <fslio.h>
#include <glib.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"

typedef struct {
    char *inputpath;
	char *maskpath;
    char *outputpath;
	char *realignmentpath;
	char *fdpath;
	char *dvarspath;
	char *flaggedpath;
	char *callString;

	double fdthreshold;
	double dvarsthreshold;

	long rpColumns;
	long rpEntries;
	double **rp;
	Point3D *maskPoints;
	unsigned long nMaskPoints;

    BOOL verbose;
    BOOL parametersValid;

    rsNiftiFile *input;
    rsNiftiFile *output;
    double ***mask;

	GOptionContext *context;

    int threads;

} rsMotionScrubbingParameters;

rsMotionScrubbingParameters *rsMotionScrubbingParseParams(int argc, char * argv[]);
rsMotionScrubbingParameters *rsMotionScrubbingInitParameters();
void rsMotionScrubbingFreeParams(rsMotionScrubbingParameters *p);
void rsMotionScrubbingPrintHelp(rsMotionScrubbingParameters *p);

#endif
