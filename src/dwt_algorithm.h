

#ifndef DWT_ALGORITHM_H
#define DWT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueDWT *DWTObj;

/***
	num radix2Exp-1 <=radix2Exp-1
	waveletType 'sym4' 'db4'/'coif4'/'fk4'/'bior4.4'
****/
int dwtObj_new(DWTObj *dwtObj,int num,int radix2Exp,
			 WaveletDiscreteType *waveletType,int *t1,int *t2);

// coefArr dataLength;mDataArr num*dataLength
void dwtObj_dwt(DWTObj dwtObj,float *dataArr,float *coefArr,float *mDataArr);

void dwtObj_free(DWTObj dwtObj);

#ifdef __cplusplus
}
#endif

#endif