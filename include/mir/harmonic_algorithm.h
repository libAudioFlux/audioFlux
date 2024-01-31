

#ifndef HARMONI_ALGORITHM_H
#define HARMONI_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueHarmonic *HarmonicObj;

/***
	samplate 32000
	radix2Exp 12 
	windowType hamm
****/
int harmonicObj_new(HarmonicObj *harmonicObj,
					int *samplate,float *lowFre,float *highFre,
					int *radix2Exp,WindowType *windowType,int *slideLength);

int harmonicObj_calTimeLength(HarmonicObj harmonicObj,int dataLength);

void harmonicObj_exec(HarmonicObj harmonicObj,float *dataArr,int dataLength);

void harmonicObj_harmonicCount(HarmonicObj harmonicObj,float low,float high,int *countArr);

void harmonicObj_free(HarmonicObj harmonicObj);

#ifdef __cplusplus
}
#endif

#endif