

#ifndef SYNSQ_ALGORITHM_H
#define SYNSQ_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueSynsq *SynsqObj;

/***
	for cwt system:cwt,dwt,wpt,st,fst
	samplate 32000
	order 1 >=1
	scaleType 
	thresh >=0 default 0.001
****/
int synsqObj_new(SynsqObj *synsqObj,int num,int radix2Exp,
				int *samplate,int *order,
				float *thresh);

void synsqObj_synsq(SynsqObj synsqObj,float *freArr,
					SpectralFilterBankScaleType scaleType,
					float *mRealArr1,float *mImageArr1,
					float *mRealArr2,float *mImageArr2);

void synsqObj_free(SynsqObj synsqObj);

#ifdef __cplusplus
}
#endif

#endif