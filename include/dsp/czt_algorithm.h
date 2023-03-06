

#ifndef CZT_ALGORITHM_H
#define CZT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueCZT *CZTObj;

int cztObj_new(CZTObj *cztObj,int radix2Exp);

void cztObj_czt(CZTObj cztObj,float *realArr1,float *imageArr1,
				float lowW,float highW,
				float *realArr3,float *imageArr3);

void cztObj_free(CZTObj cztObj);

#ifdef __cplusplus
}
#endif

#endif