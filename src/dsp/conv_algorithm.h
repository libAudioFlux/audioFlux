

#ifndef CONV_ALGORITHM_H
#define CONV_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

// convolution 相关
typedef enum{
	ConvMode_Full=0, 
	ConvMode_Same, 
	ConvMode_Valid, 

} ConvModeType;

typedef enum{
	ConvMethod_Auto=0, 
	ConvMethod_Direct, 
	ConvMethod_FFT, 

} ConvMethodType;

typedef struct OpaqueConv *ConvObj;

int convObj_new(ConvObj *convObj);

int convObj_conv(ConvObj convObj,float *vArr1,int length1,float *vArr2,int length2,
				ConvModeType *mode,ConvMethodType *method,
				float *vArr3);

void convObj_free(ConvObj convObj);


#ifdef __cplusplus
}
#endif

#endif
