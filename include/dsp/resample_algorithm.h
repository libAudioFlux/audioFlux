

#ifndef RESAMPLE_ALGORITHM_H
#define RESAMPLE_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef enum{
	ResampleAlg_Polyphase=0, // stand

	ResampleAlg_Bandlimited, // ccrma

} ResampleAlgType;

typedef enum{
	ResampleQuality_Best=0, 
	ResampleQuality_Mid, 
	ResampleQuality_Fast, 

} ResampleQualityType;


typedef struct OpaqueResample *ResampleObj;

int resampleObj_new(ResampleObj *resampleObj,ResampleQualityType *qualType,int *isScale,int *isContinue);
/***
	sinc right
	zeroNum 64, 16/32/64
	nbit 9, 5~9 echo zero-cross 1<<nbit samples
	winType hann
	value kaiser/gauss 5/2.5
	rollOff 0.945, 0.8~0.95
****/
int resampleObj_newWithWindow(ResampleObj *resampleObj,
							int *zeroNum,int *nbit,
							WindowType *winType,float *value,
							float *rollOff,
							int *isScale,
							int *isContinue);

int resampleObj_calDataLength(ResampleObj resampleObj,int dataLength);

// 32000/16000
void resampleObj_setSamplate(ResampleObj resampleObj,int sourceRate,int targetRate);
void resampleObj_setSamplateRatio(ResampleObj resampleObj,float ratio);
void resampleObj_enableContinue(ResampleObj resampleObj,int flag);

int resampleObj_resample(ResampleObj resampleObj,float *dataArr1,int dataLength1,float *dataArr2);

void resampleObj_free(ResampleObj resampleObj);
void resampleObj_debug(ResampleObj resampleObj);

#ifdef __cplusplus
}
#endif

#endif