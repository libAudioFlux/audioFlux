

#ifndef FILTERDESIGN_FREQZ_H
#define FILTERDESIGN_FREQZ_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
/***
	滤波器设计相关模块
	1. FIR/IIR
	2. 经典方法
	3. 特殊场景使用 gammatone/mel/bark等
****/

void filterDesign_freqzBA(float *bArr,float *aArr,int length,
						int fftLength,int samplate,int isWhole,
						float *kArr,
						float *realArr,float *imageArr,
						float *wArr);

void filterDesign_freqzSOS(float *mArr,int nLength,
						int fftLength,int samplate,int isWhole,
						float *kArr,
						float *realArr,float *imageArr,
						float *wArr);

/***
	b[0]+b[1]*e^(-jw)+...+b[M]*e^(-j(M-1)w);
	wArr1 0~2*pi 刻度频率
****/
void filterDesign_calFreResponse(float *wArr1,int length1,
							float *vArr2,int length2,
							float *realArr3,float *imageArr3);



#ifdef __cplusplus
}
#endif

#endif