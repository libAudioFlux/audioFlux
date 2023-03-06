

#ifndef HHT_ALGORITHM_H
#define HHT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueHHT *HHTObj;

/***
	num radix2Exp-1 <=radix2Exp-1
	waveletType 'sym4' 'db4'/'coif4'/'fk4'/'bior4.4'
****/
int hhtObj_new(HHTObj *hhtObj,int num,int radix2Exp,
			 WaveletDiscreteType *waveletType,int *t1,int *t2);

// coefArr dataLength;mDataArr num*dataLength
void hhtObj_hht(HHTObj hhtObj,float *dataArr,float *coefArr,float *mDataArr);

void hhtObj_free(HHTObj hhtObj);

#ifdef __cplusplus
}
#endif

#endif