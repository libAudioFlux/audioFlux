

#ifndef DFT_ALGORITHM_H
#define DFT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueDFT *DFTObj;

int dftObj_new(DFTObj *dftObj,int length);

void dftObj_dft(DFTObj dftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2);
void dftObj_idft(DFTObj dftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2);

void dftObj_free(DFTObj dftObj);

#ifdef __cplusplus
}
#endif

#endif