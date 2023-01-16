

#ifndef HARMONICRATIO_ALGORITHM_H
#define HARMONICRATIO_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueHarmonicRatio *HarmonicRatioObj;

/***
	samplate 32000
	lowFre 25
	radix2Exp 12 
	windowType hamm
****/
int harmonicRatioObj_new(HarmonicRatioObj *harmonicRatioObj,
						int *samplate,float *lowFre,
						int *radix2Exp,WindowType *windowType,int *slideLength);

int harmonicRatioObj_calTimeLength(HarmonicRatioObj harmonicRatioObj,int dataLength);
void harmonicRatioObj_harmonicRatio(HarmonicRatioObj harmonicRatioObj,float *dataArr,int dataLength,float *valueArr);

void harmonicRatioObj_free(HarmonicRatioObj harmonicRatioObj);

#ifdef __cplusplus
}
#endif

#endif