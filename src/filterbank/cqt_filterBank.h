

#ifndef CQT_FILTERBANK_H
#define CQT_FILTERBANK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

// filterBank相关
void cqt_filterBank(int num,float *freBandArr,int samplate,
				int binPerOctave,SpectralFilterBankNormalType normType,WindowType *winType,
				float *factor,float *beta,float *thresh,
				float *lenArr,int fftLength,
				float *mRealFilterBankArr,float *mImageFilterBankArr);

void cqt_downFilterBank(int num,float *freBandArr,int samplate,
				int binPerOctave,SpectralFilterBankNormalType normType,WindowType *winType,
				float *factor,float *beta,float *thresh,
				float *lenArr,int fftLength,
				float *mRealFilterBankArr,float *mImageFilterBankArr);

float cqt_calQ(int binPerOctave,float factor);
float *cqt_calFreArr(float minFre,int num,int binPerOctave);

int cqt_calFFTLength(float minFre,int samplate,
				int binPerOctave,
				float *factor,float *beta);

void cqt_calLengthArr(int num,float *freBandArr,int samplate,
				int binPerOctave,
				float *factor,float *beta,
				float *lenArr);

#ifdef __cplusplus
}
#endif

#endif