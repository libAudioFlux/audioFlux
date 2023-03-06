

#ifndef EMD_ALGORITHM_H
#define EMD_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueEMD *EMDObj;

int emdObj_new(EMDObj *emdObj,int num,int radix2Exp);

void emdObj_emd(EMDObj emdObj,float *dataArr,float *mRealArr1,float *mImageArr1);

void emdObj_free(EMDObj emdObj);

#ifdef __cplusplus
}
#endif

#endif