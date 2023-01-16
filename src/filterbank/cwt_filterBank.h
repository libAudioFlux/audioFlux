// clang 

#ifndef CWT_FILTERBANK_H
#define CWT_FILTERBANK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

/***
	cwt-filterBank相关 

	log+linear+mel/bark/erb 都不包含边界
		log 针对cwt-spectrogram形式上包含
		log is use default

	Morse gamma=3 beta=20; 
	Morlet gamma=6 beta=2; 
	Bump gamma=5 beta=0.6; gamma(3,6),beta(0.1,1.2)
	
	Paul gamma=4
	DOG gamma=2 beta=2;
	Mexican beta=2

	type semantic is wavelet type and window type(styleType)
****/
void cwt_filterBank(int num,int dataLength,int samplate,int padLength,
					WaveletContinueType type,float gamma,float beta,
					SpectralFilterBankScaleType scaleType,
					float lowFre,float highFre,int binPerOctave,
					float *mFilterBankArr,
					float *freBandArr,
					int *binBandArr);



#ifdef __cplusplus
}
#endif

#endif