

#ifndef CEPSTROGRAM_ALGORITHM_H
#define CEPSTROGRAM_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueCepstrogram *CepstrogramObj;

/***
	radix2Exp 12 
	WindowType "rect"
	slideLength 1024
****/
int cepstrogramObj_new(CepstrogramObj *cepstrogramObj,int radix2Exp,WindowType *windowType,int *slideLength);

int cepstrogramObj_calTimeLength(CepstrogramObj cepstrogramObj,int dataLength);

/***
	cepNum 4~128 ,formant estimate number
	mDataArr1 cepstrums ,timeLength*(fftLength/2+1)
	mDataArr2 envelope(formant) ,timeLength*(fftLength/2+1)
	mDataArr3 details(tone) ,timeLength*(fftLength/2+1)
****/
void cepstrogramObj_cepstrogram(CepstrogramObj cepstrogramObj,int cepNum,float *dataArr,int dataLength,
								float *mDataArr1,float *mDataArr2,float *mDataArr3);

void cepstrogramObj_cepstrogram2(CepstrogramObj cepstrogramObj,int cepNum,float *mRealArr,float *mImageArr,int nLength,
								float *mDataArr1,float *mDataArr2,float *mDataArr3);

void cepstrogramObj_enableDebug(CepstrogramObj cepstrogramObj,int flag);

void cepstrogramObj_free(CepstrogramObj cepstrogramObj);

#ifdef __cplusplus
}
#endif

#endif