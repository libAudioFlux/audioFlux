

#ifndef BFT_ALGORITHM_H
#define BFT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueBFT *BFTObj;

/***
	num>=2&&num<=2048
	samplate 32000
	binPerOctave 12 >=4&&<=48
	
	radix2Exp 12 
	WindowType "hann"
	slideLength 1024

	filterScaleType "linear"
	filterStyleType "slaney"
	filterNormalType "none"
	dataType "power"

	isReassign 0
	isTemporal 0
****/
int bftObj_new(BFTObj *bftObj,int num,int radix2Exp,
			int *samplate,float *lowFre,float *highFre,int *binPerOctave,
			WindowType *windowType,int *slideLength,
			SpectralFilterBankScaleType *filterScaleType,
			SpectralFilterBankStyleType *filterStyleType,
			SpectralFilterBankNormalType *filterNormalType,
			SpectralDataType *dataType,
			int *isReassign,
			int *isTemporal);

int bftObj_calTimeLength(BFTObj bftObj,int dataLength);

float *bftObj_getFreBandArr(BFTObj bftObj);
int *bftObj_getBinBandArr(BFTObj bftObj);

// 0 complex,1 real ->mRealArr
void bftObj_setResultType(BFTObj bftObj,int type);
void bftObj_setDataNormValue(BFTObj bftObj,float normValue);

void bftObj_bft(BFTObj bftObj,float *dataArr,int dataLength,float *mRealArr3,float *mImageArr3);

// energy/rms/zeroCrossRate
void bftObj_getTemporalData(BFTObj bftObj,float **eArr,float **rArr,float **zArr);

void bftObj_free(BFTObj bftObj);

#ifdef __cplusplus
}
#endif

#endif