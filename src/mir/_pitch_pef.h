

#ifndef _PITCH_PEF_H
#define _PITCH_PEF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaquePitchPEF *PitchPEFObj;

/***
	samplate 32000
	lowFre 32,
	highFre 2000
	cutFre 4000, >highFre

	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	alpha 10 >0, beta 0.5 0~1, gamma 1.8 >1
	
	isContinue 0
****/
int pitchPEFObj_new(PitchPEFObj *pitchPEFObj,
				int *samplate,float *lowFre,float *highFre,float *cutFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				float *alpha,float *beta,float *gamma,
				int *isContinue);

int pitchPEFObj_calTimeLength(PitchPEFObj pitchPEFObj,int dataLength);
void pitchPEFObj_setFilterParams(PitchPEFObj pitchPEFObj,float alpha,float beta,float gamma);

void pitchPEFObj_pitch(PitchPEFObj pitchPEFObj,float *dataArr,int dataLength,
					float *freArr);

void pitchPEFObj_enableDebug(PitchPEFObj pitchPEFObj,int isDebug);
void pitchPEFObj_free(PitchPEFObj pitchPEFObj);

#ifdef __cplusplus
}
#endif

#endif