

#ifndef FILTERDESIGN_FIR_H
#define FILTERDESIGN_FIR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

// 窗函数设计法 默认hamm
float *filterDesign_fir1(int order,
					float *wcArr,
					FilterBandType bandType,
					WindowType *winType,float *value,
					int *isNoScale);

float *filterDesign_fir2(int order,
					float *wcArr,
					FilterBandType bandType,
					float *winArr,
					int *isNoScale);

// order must odd; first derivative
float *filterDesign_smooth1(int order);
float *filterDesign_mean(int order);

// filter(b,a,x)
void filterDesign_filter(float *bArr,float *aArr,float *xArr,
						int bLength,int aLength,int xLength,
						float *yArr);

// filtfilt(b,a,x)
void filterDesign_filtfilt(float *bArr,float *aArr,float *xArr,
						int bLength,int aLength,int xLength,
						float *yArr);



#ifdef __cplusplus
}
#endif

#endif