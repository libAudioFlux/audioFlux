

#ifndef NSGT_FILTERBANK_H
#define NSGT_FILTERBANK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

/***
	nsgt-filterBank相关 
	log+linear+mel/bark/erb 都不包含边界
		log 针对nsgt-spectrogram形式上包含
	log is nsg-cqt
	type 0 efficient 1 standard 同论文
****/
void nsgt_filterBank(int num,int dataLength,int samplate, int minLength,
					int type,
					SpectralFilterBankScaleType scaleType,
					SpectralFilterBankStyleType styleType,
					SpectralFilterBankNormalType normType,
					float lowFre,float highFre,int binPerOctave,
					float **filterBankArr,int *lengthArr,
					float *freBandArr,int *binBandArr,int *offsetArr,
					int *maxLength,int *totalLength);


#ifdef __cplusplus
}
#endif

#endif
