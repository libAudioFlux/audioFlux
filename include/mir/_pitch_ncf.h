

#ifndef _PITCH_NCF_H
#define _PITCH_NCF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaquePitchNCF *PitchNCFObj;

/***
	samplate 32000
	lowFre 32,
	highFre 2000
	radix2Exp 12
	WindowType rect
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int pitchNCFObj_new(PitchNCFObj *pitchNCFObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				int *isContinue);

int pitchNCFObj_calTimeLength(PitchNCFObj pitchNCFObj,int dataLength);

void pitchNCFObj_pitch(PitchNCFObj pitchNCFObj,float *dataArr,int dataLength,
					float *freArr);

void pitchNCFObj_enableDebug(PitchNCFObj pitchNCFObj,int isDebug);
void pitchNCFObj_free(PitchNCFObj pitchNCFObj);

#ifdef __cplusplus
}
#endif

#endif