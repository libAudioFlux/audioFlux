

#ifndef _PITCH_CEP_H
#define _PITCH_CEP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaquePitchCEP *PitchCEPObj;

/***
	samplate 32000
	lowFre 32,
	highFre 2000
	radix2Exp 12
	WindowType Hamm
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int pitchCEPObj_new(PitchCEPObj *pitchCEPObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				int *isContinue);

int pitchCEPObj_calTimeLength(PitchCEPObj pitchCEPObj,int dataLength);

void pitchCEPObj_pitch(PitchCEPObj pitchCEPObj,float *dataArr,int dataLength,
					float *freArr);

void pitchCEPObj_enableDebug(PitchCEPObj pitchCEPObj,int isDebug);
void pitchCEPObj_free(PitchCEPObj pitchCEPObj);

#ifdef __cplusplus
}
#endif

#endif