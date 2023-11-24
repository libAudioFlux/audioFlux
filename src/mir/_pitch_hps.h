

#ifndef _PITCH_HPS_H
#define _PITCH_HPS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaquePitchHPS *PitchHPSObj;

/***
	samplate 32000
	lowFre 32,
	highFre 2000
	radix2Exp 12
	WindowType Hamm
	slideLength (1<<radix2Exp)/4
	
	harmonicCount 5 >0
	isContinue 0
****/
int pitchHPSObj_new(PitchHPSObj *pitchHPSObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				int *harmonicCount,
				int *isContinue);

int pitchHPSObj_calTimeLength(PitchHPSObj pitchHPSObj,int dataLength);

void pitchHPSObj_pitch(PitchHPSObj pitchHPSObj,float *dataArr,int dataLength,
					float *freArr);

void pitchHPSObj_enableDebug(PitchHPSObj pitchHPSObj,int isDebug);
void pitchHPSObj_free(PitchHPSObj pitchHPSObj);

#ifdef __cplusplus
}
#endif

#endif