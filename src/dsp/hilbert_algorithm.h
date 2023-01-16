

#ifndef HILBERT_ALGORITHM_H
#define HILBERT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueHilbert *HilbertObj;

int hilbertObj_new(HilbertObj *hilbertObj,int radix2Exp);

void hilbertObj_hilbert(HilbertObj hilbertObj,float *dataArr,
						float *realArr3,float *imageArr3);

void hilbertObj_free(HilbertObj hilbertObj);

#ifdef __cplusplus
}
#endif

#endif