

#ifndef DECONV_ALGORITHM_H
#define DECONV_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueDeconv *DeconvObj;

int deconvObj_new(DeconvObj *deconvObj,int num);

void deconvObj_setTimeLength(DeconvObj deconvObj,int timeLength);

/***
	mDataArr1 'mag' mag||power
	mDataArr2 timbre(formant) 
	mDataArr3 pitch
****/
void deconvObj_deconv(DeconvObj deconvObj,float *mDataArr1,float *mRealArr3,float *mImageArr3);

void deconvObj_free(DeconvObj deconvObj);

#ifdef __cplusplus
}
#endif

#endif