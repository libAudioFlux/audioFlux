

#ifndef TEMPORAL_ALGORITHM_H
#define TEMPORAL_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueTemporal *TemporalObj;

/***
	frameLength 2048
	slideLength 512
	windowType Hann
****/
int temporalObj_new(TemporalObj *temporalObj,
			 		int *frameLength,int *slideLength,WindowType *windowType);

int temporalObj_calTimeLength(TemporalObj temporalObj,int dataLength);

void temporalObj_temporal(TemporalObj temporalObj,float *dataArr,int dataLength);

/***
	energy/rms/zeroCrossRate
	mDataArr timeLength*frameLength
****/
void temporalObj_getData(TemporalObj temporalObj,
						float **eArr,float **rArr,float **zArr,
						float **mDataArr);

// gamma 1/10/20 song 0.5
void temporalObj_ezr(TemporalObj temporalObj,float gamma,float *vArr3);

void temporal_free(TemporalObj temporalObj);

#ifdef __cplusplus
}
#endif

#endif