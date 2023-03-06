

#ifndef FST_ALGORITHM_H
#define FST_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueFST *FSTObj;

int fstObj_new(FSTObj *fstObj,int radix2Exp);

void fstObj_fst(FSTObj fstObj,float *dataArr,int minIndex,int maxIndex,float *mRealArr,float *mImageArr);

void fstObj_free(FSTObj fstObj);


#ifdef __cplusplus
}
#endif

#endif