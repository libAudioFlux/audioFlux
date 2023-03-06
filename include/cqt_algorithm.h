

#ifndef CQT_ALGORITHM_H
#define CQT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueCQT *CQTObj;

int cqtObj_new(CQTObj *cqtObj,int num,int samplate,float minFre,int *isContinue);
/***
	num n*binPerOctave
	samplate 32000
	minFre C1=32.703
	binPerOctave 12
	factor 1 
		>0  <1 short window
	beta 0 
		>=0 0 cqt >0 vqt ???
	thresh 0.01
	windowType 'hann'
	slideLength fftLength/4
	isScale 1
****/
int cqtObj_newWith(CQTObj *cqtObj,int num,
				int *samplate,float *minFre,int *binPerOctave,
				float *factor,float *beta,float *thresh,
				WindowType *windowType,int *slideLength,int *isContinue,
				SpectralFilterBankNormalType *filterNormalType,int *isScale);

int cqtObj_calTimeLength(CQTObj cqtObj,int dataLength);
int cqtObj_getFFTLength(CQTObj cqtObj);
// num
float *cqtObj_getFreBandArr(CQTObj cqtObj);

void cqtObj_setScale(CQTObj cqtObj,int flag);

void cqtObj_cqt(CQTObj cqtObj,float *dataArr,int dataLength,float *mRealArr,float *mImageArr);
/***
	chromaNum 12
	dataType 'power' mag||power 
	mDataArr timeLength*chromaNum
	normType 'Max'
****/
void cqtObj_chroma(CQTObj cqtObj,int *chromaNum,SpectralDataType *dataType,ChromaDataNormalType *normType,
				float *mRealArr,float *mImageArr,
				float *mDataArr);

// mDataArr1  "mag" mag||power
void cqtObj_cqcc(CQTObj cqtObj,float *mDataArr1,int ccNum,CepstralRectifyType *rectifyType,float *mDataArr2);

// mDataArr1  "power" mag||power
void cqtObj_cqhc(CQTObj cqtObj,float *mDataArr1,int hcNum,float *mDataArr2);
// mDataArr1  "mag" mag||power
void cqtObj_deconv(CQTObj cqtObj,float *mDataArr1,float *mDataArr2,float *mDataArr3);

void cqtObj_free(CQTObj cqtObj);



#ifdef __cplusplus
}
#endif

#endif