

#ifndef CHROMA_FILTERBANK_H
#define CHROMA_FILTERBANK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

// filterBank相关
/***
	octaveCenter 5 A440
	octaveWidth 2 
****/
void chroma_stftFilterBank(int num,int fftLength,int samplate,
						float *octaveCenter,float *octaveWidth,
						float *mFilterBankArr);

// minFre C1=32.703
void chroma_cqtFilterBank(int num,int cqtLength,int binPerOctave,
						float *minFre,
						float *mFilterBankArr);

void chroma_genericFilterBank(int num,int fftLength,int samplate,
							int freLength,float *freBandArr,
							float *mFilterBankArr);

#ifdef __cplusplus
}
#endif

#endif