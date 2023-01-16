

#ifndef STFT_ALGORITHM_H
#define STFT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueSTFT *STFTObj;

int stftObj_new(STFTObj *stftObj,int radix2Exp,WindowType *windowType,int *slideLength,int *isContinue);

void stftObj_setSlideLength(STFTObj stftObj,int slideLength);
// value1 => constant/left
void stftObj_setPadding(STFTObj stftObj,
					PaddingPositionType *positionType,PaddingModeType *modeType,
					float *value1,float *value2);

void stftObj_useWindowDataArr(STFTObj stftObj,float *winDataArr);
float *stftObj_getWindowDataArr(STFTObj stftObj);

// center zero padding
void stftObj_enablePadding(STFTObj stftObj,int flag);
void stftObj_enableContinue(STFTObj stftObj,int flag);

// set/enable后执行
int stftObj_calTimeLength(STFTObj stftObj,int dataLength);
int stftObj_calDataLength(STFTObj stftObj,int timeLength);

void stftObj_stft(STFTObj stftObj,float *dataArr,int dataLength,float *mRealArr,float *mImageArr);
// type 0(deault) 'weight overlap-add' 1 'overlap-add'
void stftObj_istft(STFTObj stftObj,float *mRealArr,float *mImageArr,int nLength,int type,float *dataArr);

void stftObj_free(STFTObj stftObj);
void stftObj_debug(STFTObj stftObj);

#ifdef __cplusplus
}
#endif

#endif