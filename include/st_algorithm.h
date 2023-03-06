

#ifndef ST_ALGORITHM_H
#define ST_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueST *STObj;

// minIndex>0 maxIndex<fftLength/2 not contain edge
int stObj_new(STObj *stObj,int radix2Exp,int minIndex,int maxIndex,float *factor,float *norm);

void stObj_useBinArr(STObj stObj,int *binArr,int length);
void stObj_setValue(STObj stObj,float factor,float norm);

void stObj_st(STObj stObj,float *dataArr,float *mRealArr,float *mImageArr);

void stObj_free(STObj stObj);


#ifdef __cplusplus
}
#endif

#endif