#ifndef rstools_rscorrelation_common_h
#define rstools_rscorrelation_common_h

#include "rscorrelation_ui.h"

#ifdef __cplusplus
extern "C" {
#endif

void rsCorrelationInit(rsCorrelationParameters* p);
void rsCorrelationRun(rsCorrelationParameters *p);
void rsCorrelationDestroy(rsCorrelationParameters* p);

void rsCorrelationWriteCorrelationFile(rsCorrelationParameters* p);
double *rsReadRegressorFromStream(FILE *stream, unsigned int *nValues);

#ifdef __cplusplus
}
#endif

#endif