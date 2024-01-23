

#ifndef PHASE_VOCODER_H
#define PHASE_VOCODER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

/***
	mDataArr1 stft time*fftLength
	rate 0.5~2
	slideLength fftLength/4
****/
void phase_vocoder(float *mRealArr1,float *mImageArr1,int timeLength,int fftLength,int slideLength,float rate,
				float *mRealArr2,float *mImageArr2);


#ifdef __cplusplus
}
#endif

#endif