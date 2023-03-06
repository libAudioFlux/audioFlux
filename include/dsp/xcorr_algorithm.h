

#ifndef XCORR_ALGORITHM_H
#define XCORR_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef enum{
	XcorrNormal_None=0, 
	XcorrNormal_Coeff, 

} XcorrNormalType;

typedef struct OpaqueXcorr *XcorrObj;

int xcorrObj_new(XcorrObj *xcorrObj);

/***
	vArr1 vArr2 xcorr 
	vArr1 NULL	autocorr 
		vArr1==vArr2
	normType default 'Coeff'
	return
		vArr3 2*length-1
		index vArr3 actual is index-(length-1)
****/
int xcorrObj_xcorr(XcorrObj xcorrObj,float *vArr1,float *vArr2,int length,
				XcorrNormalType *normType,
				float *vArr3,float *maxValue);

void xcorrObj_free(XcorrObj xcorrObj);


#ifdef __cplusplus
}
#endif

#endif
