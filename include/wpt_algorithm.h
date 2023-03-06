

#ifndef WPT_ALGORITHM_H
#define WPT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueWPT *WPTObj;

/***
	num radix2Exp-1 <=radix2Exp-1
	waveletType 'sym4' 'db4'/'coif4'/'fk4'/'bior4.4'
****/
int wptObj_new(WPTObj *wptObj,int num,int radix2Exp,
			 WaveletDiscreteType *waveletType,int *t1,int *t2);

// coefArr dataLength;mDataArr 2^num*dataLength
void wptObj_wpt(WPTObj wptObj,float *dataArr,float *coefArr,float *mDataArr);

void wptObj_free(WPTObj wptObj);

#ifdef __cplusplus
}
#endif

#endif