

#ifndef _PITCH_LHS_H
#define _PITCH_LHS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaquePitchLHS *PitchLHSObj;

/***
	samplate 32000
	lowFre 32,
	highFre 2000

	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	
	harmonicCount 5 >0
	isContinue 0
****/
int pitchLHSObj_new(PitchLHSObj *pitchLHSObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				int *harmonicCount,
				int *isContinue);

int pitchLHSObj_calTimeLength(PitchLHSObj pitchLHSObj,int dataLength);

void pitchLHSObj_pitch(PitchLHSObj pitchLHSObj,float *dataArr,int dataLength,
					float *freArr);

void pitchLHSObj_enableDebug(PitchLHSObj pitchLHSObj,int isDebug);
void pitchLHSObj_free(PitchLHSObj pitchLHSObj);

#ifdef __cplusplus
}
#endif

#endif