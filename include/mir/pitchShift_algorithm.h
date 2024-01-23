

#ifndef PITCHSHIFT_ALGORITHM_H
#define PITCHSHIFT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaquePitchShift *PitchShiftObj;

/***
	radix2Exp 12
	WindowType hann
	slideLength (1<<radix2Exp)/4
****/
int pitchShiftObj_new(PitchShiftObj *pitchShiftObj,int *radix2Exp,int *slideLength,WindowType *windowType);

// nSemitone -12~12
void pitchShiftObj_pitchShift(PitchShiftObj pitchShiftObj,int samplate,int nSemitone,float *dataArr1,int dataLength1,float *dataArr2);

void pitchShiftObj__free(PitchShiftObj pitchShiftObj);

#ifdef __cplusplus
}
#endif

#endif