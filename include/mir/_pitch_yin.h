

#ifndef _PITCH_YIN_H
#define _PITCH_YIN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

typedef struct OpaquePitchYIN *PitchYINObj;

/***
	samplate 32000
	lowFre 27
	highFre 2000

	radix2Exp 12
	slideLength (1<<radix2Exp)/4
	autoLength (1<<radix2Exp)/2

	isContinue 0
****/
int pitchYINObj_new(PitchYINObj *pitchYINObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,int *autoLength,
				int *isContinue);

// default 0.6 thresh>0&&thresh<1
void pitchYINObj_setThresh(PitchYINObj pitchYINObj,float thresh);
int pitchYINObj_calTimeLength(PitchYINObj pitchYINObj,int dataLength);

void pitchYINObj_pitch(PitchYINObj pitchYINObj,float *dataArr,int dataLength,
					float *freArr,float *valueArr1,float *valueArr2);

int pitchYINObj_getTroughData(PitchYINObj pitchYINObj,float **mFreArr,float **mTroughArr,int **lenArr);

void pitchYINObj_enableDebug(PitchYINObj pitchYINObj,int isDebug);
void pitchYINObj_free(PitchYINObj pitchYINObj);

#ifdef __cplusplus
}
#endif

#endif