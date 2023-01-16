

#ifndef WSST_ALGORITHM_H
#define WSST_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueWSST *WSSTObj;

/***
	waveletType 'morlet'
	scaleType 'log'

	'morse' gamma 3 beta 20
	'morlet' gamma 6 beta 2
	'bump' gamma 5 beta 0.6

	'paul' gamma 4
	'dog' gamma 2 beta 2; must even
	'mexican' beta 2 

	thresh >=0 default 0.001
****/
int wsstObj_new(WSSTObj *wsstObj,int num,int radix2Exp,
				int *samplate,float *lowFre,float *highFre,int *binPerOctave,
			 	WaveletContinueType *waveletType,
			 	SpectralFilterBankScaleType *scaleType,
			 	float *gamma,float *beta,
			 	float *thresh,
			 	int *isPadding);

float *wsstObj_getFreBandArr(WSSTObj wsstObj);
int *wsstObj_getBinBandArr(WSSTObj wsstObj);

// order >=1
void wsstObj_setOrder(WSSTObj wsstObj,int order);

// mRealArr2&mImageArr2 may NULL
void wsstObj_wsst(WSSTObj wsstObj,float *dataArr,
				float *mRealArr1,float *mImageArr1,
				float *mRealArr2,float *mImageArr2);

void wsstObj_free(WSSTObj wsstObj);

#ifdef __cplusplus
}
#endif

#endif