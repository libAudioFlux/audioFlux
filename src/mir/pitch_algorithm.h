

#ifndef PITCH_ALGORITHM_H
#define PITCH_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

typedef enum{
	Pitch_YIN=0, 
	Pitch_STFT,
	Pitch_NCF,
	Pitch_PEF,

} PitchType;

typedef struct OpaquePitch *PitchObj;

/***
	type None
	samplate 32000
	lowFre 27
	highFre 2000

	radix2Exp 12
	slideLength (1<<radix2Exp)/4
	autoLength (1<<radix2Exp)/2

	isContinue 0
****/
int pitchObj_new(PitchObj *pitchObj,PitchType *type,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,int *autoLength,
				int *isContinue);

// default 0.1 thresh>0&&thresh<1
void pitchObj_setThresh(PitchObj pitchObj,float thresh);
int pitchObj_calTimeLength(PitchObj pitchObj,int dataLength);

void pitchObj_pitch(PitchObj pitchObj,float *dataArr,int dataLength,
					float *freArr,float *valueArr1,float *valueArr2);

void pitchObj_enableDebug(PitchObj pitchObj,int isDebug);
void pitchObj_free(PitchObj pitchObj);

#ifdef __cplusplus
}
#endif

#endif