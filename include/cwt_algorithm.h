

#ifndef CWT_ALGORITHM_H
#define CWT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueCWT *CWTObj;

/***
	waveletType 'morse'
	scaleType 'log'
	'morse' gamma 3 beta 20
	'morlet' gamma 6 beta 2
	'bump' gamma 5 beta 0.6

	'paul' gamma 4
	'dog' gamma 2 beta 2; must even
	'mexican' beta 2 
	
	'hermit' gamma 5 beta 2
	'ricker' gamma 4
****/
int cwtObj_new(CWTObj *cwtObj,int num,int radix2Exp,
			 int *samplate,float *lowFre,float *highFre,int *binPerOctave,
			 WaveletContinueType *waveletType,
			 SpectralFilterBankScaleType *scaleType,
			 float *gamma,float *beta,
			 int *isPadding);

float *cwtObj_getFreBandArr(CWTObj cwtObj);
int *cwtObj_getBinBandArr(CWTObj cwtObj);

void cwtObj_cwt(CWTObj cwtObj,float *dataArr,float *mRealArr3,float *mImageArr3);

void cwtObj_enableDet(CWTObj cwtObj,int flag);
void cwtObj_cwtDet(CWTObj cwtObj,float *dataArr,float *mRealArr3,float *mImageArr3);

void cwtObj_free(CWTObj cwtObj);

#ifdef __cplusplus
}
#endif

#endif