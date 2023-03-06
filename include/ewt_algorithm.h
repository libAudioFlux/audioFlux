

#ifndef EWT_ALGORITHM_H
#define EWT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueEWT *EWTObj;

/***
	num radix2Exp-1 <=radix2Exp-1
	waveletType 'sym4' 'db4'/'coif4'/'fk4'/'bior4.4'
****/
int ewtObj_new(EWTObj *ewtObj,int num,int radix2Exp,
			 WaveletDiscreteType *waveletType,int *t1,int *t2);

// coefArr dataLength;mDataArr num*dataLength
void ewtObj_ewt(EWTObj ewtObj,float *dataArr,float *coefArr,float *mDataArr);

void ewtObj_free(EWTObj ewtObj);

#ifdef __cplusplus
}
#endif

#endif