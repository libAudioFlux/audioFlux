

#ifndef XXCC_ALGORITHM_H
#define XXCC_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueXXCC *XXCCObj;

int xxccObj_new(XXCCObj *xxccObj,int num);

void xxccObj_setTimeLength(XXCCObj xxccObj,int timeLength);

/***
	mDataArr1 'mag' mag||power
	ccNum<=num
	rectifyType CepstralRectify_Log
****/
void xxccObj_xxcc(XXCCObj xxccObj,float *mDataArr1,int ccNum,CepstralRectifyType *rectifyType,float *mDataArr2);

/***
	mfcc standard/xxcc standard
	mDataArr1 'mag' mag||power
	ccNum<=num
	deltaWindowLength 9(defalut); must odd>=3
	energyType Repalce; Append ccNum+1
	rectifyType Log
	return timeLength*(ccNum+1?)
****/
void xxccObj_xxccStandard(XXCCObj xxccObj,float *mDataArr1,int ccNum,float *energyArr,
						int *deltaWindowLength,CepstralEnergyType *energyType,CepstralRectifyType *rectifyType,
						float *mCoeArr,float *mDeltaArr1,float *mDeltaArr2);

void xxccObj_free(XXCCObj xxccObj);

#ifdef __cplusplus
}
#endif

#endif