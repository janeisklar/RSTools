#ifndef rstools_rscorrelation_common_h
#define rstools_rscorrelation_common_h

#include "rscorrelation_ui.h"

void rsCorrelationInit(rsCorrelationParameters* p);
void rsCorrelationRun(rsCorrelationParameters *p);
void rsCorrelationDestroy(rsCorrelationParameters* p);

void rsCorrelationWriteCorrelationFile(rsCorrelationParameters* p);
double* rsReadRegressorFromStandardInput(unsigned int *nValues);

#endif