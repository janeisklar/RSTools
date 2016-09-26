#ifndef rstools_scrubbing_common_h
#define rstools_scrubbing_common_h

#include "rssmoothing_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsSmoothingInit(rsSmoothingParameters *p);
void rsSmoothingRun(rsSmoothingParameters *p);
void rsSmoothingDestroy(rsSmoothingParameters *p);

double ***rsCreateGaussianKernel(double sigma, short *xdimkernel, short *ydimkernel, short *zdimkernel, double xvoxsize, double yvoxsize, double zvoxsize);
void rsConvolveWithKernel(double ***result, double ***input, double ***kernel, short xdim, short ydim, short zdim, short xdimKernel, short ydimKernel, short zdimKernel);

#ifdef __cplusplus
}
#endif

#endif
