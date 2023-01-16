

#ifndef _PITCH_STFT_H
#define _PITCH_STFT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaquePitchSTFT *PitchSTFTObj;

/***
	samplate 32000
	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int pitchSTFTObj_new(PitchSTFTObj *pitchSTFTObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				int *isContinue);

int pitchSTFTObj_calTimeLength(PitchSTFTObj pitchSTFTObj,int dataLength);

void pitchSTFTObj_pitch(PitchSTFTObj pitchSTFTObj,float *dataArr,int dataLength,
					float *freArr,float *dbArr);

int pitchSTFTObj_getCorrData(PitchSTFTObj pitchSTFTObj,float **mCorr,int **lenArr);

void pitchSTFTObj_enableDebug(PitchSTFTObj pitchSTFTObj,int isDebug);
void pitchSTFTObj_free(PitchSTFTObj pitchSTFTObj);

#ifdef __cplusplus
}
#endif

#endif