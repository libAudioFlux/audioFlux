

#ifndef SWT_ALGORITHM_H
#define SWT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueSWT *SWTObj;

/***
	num >=1
	fftLength fftLength%2^num=0
	waveletType 'sym4' 'db4'/'coif4'/'fk4'/'bior4.4'
****/
int swtObj_new(SWTObj *swtObj,int num,int fftLength,
			 WaveletDiscreteType *waveletType,int *t1,int *t2);

// num*dataLength,mDataArr1 app,mDataArr2 det
void swtObj_swt(SWTObj swtObj,float *dataArr,float *mDataArr1,float *mDataArr2);

void swtObj_free(SWTObj swtObj);

#ifdef __cplusplus
}
#endif

#endif