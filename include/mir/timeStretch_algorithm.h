

#ifndef TIMESTRETCH_ALGORITHM_H
#define TIMESTRETCH_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaqueTimeStretch *TimeStretchObj;

/***
	radix2Exp 12
	WindowType hann
	slideLength (1<<radix2Exp)/4
****/
int timeStretchObj_new(TimeStretchObj *timeStretchObj,int *radix2Exp,int *slideLength,WindowType *windowType);

int timeStretchObj_calDataCapacity(TimeStretchObj timeStretchObj,float rate,int dataLength);
int timeStretchObj_timeStretch(TimeStretchObj timeStretchObj,float rate,float *dataArr1,int dataLength1,float *dataArr2);

void timeStretchObj_free(TimeStretchObj timeStretchObj);

#ifdef __cplusplus
}
#endif

#endif