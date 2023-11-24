

#ifndef FLUX_UTIL_H
#define FLUX_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

// 2^n相关
int util_isPowerTwo(int value);
int util_ceilPowerTwo(int value);
int util_floorPowerTwo(int value);
int util_roundPowerTwo(int value);

// value=2^n =>n
int util_powerTwoBit(int value);

// 最大公约数 a>b
int util_gcd(int a,int b);

// fre相关
float util_midiToFre(int midi);
int util_freToMidi(float fre);
int util_midiTimes(int midi1,int midi2);
int util_freTimes(float fre1,float fre2);

int util_freToSimularMidi(float fre);
int util_freTimes1(float fre1,float fre2);

// scale
void util_minMaxScale(float *vArr1,int length,float *vArr2);
void util_standScale(float *vArr1,int length,int type,float *vArr2);
void util_maxAbsScale(float *vArr1,int length,float *vArr2);
void util_robustScale(float *vArr1,int length,float *vArr2);

void util_centerScale(float *vArr1,int length,float *vArr2);
void util_meanScale(float *vArr1,int length,float *vArr2);
void util_arctanScale(float *vArr1,int length,float *vArr2);

// 向量归一化 type 0 p 1/2/3... 1 Inf 2 -Inf
void util_normalize(float *vArr1,int length,int type,float p,float *vArr2);

// mag/power/DB min=-80 
void util_powerToDB(float *pArr,int length,float min,float *dArr);

void util_powerToAbsDB(float *pArr,int length,int fftLength,int isNorm,float min,float *dArr);
void util_magToAbsDB(float *pArr,int length,int fftLength,int isNorm,float min,float *dArr);

// log compression gamma=1.0 1/10/20/...
void util_logCompress(float *vArr1,float *gamma,int length,float *vArr2);
void util_log10Compress(float *vArr1,float *gamma,int length,float *vArr2);

// gamma eps=1e-5
float util_gamma(float x,float *eps);
long double util_gammal(long double x,long double *eps);

// quadratically interpolated return p ∈[-1/2,1/2]
float util_qaudInterp(float value1,float value2,float value3,float *outValue);

// peak pick
void util_peakPick(float *pArr,int length,int start,int end,int distance,int num,float *vArr,int *cArr);

// order must odd; delta/deltaDelta
void util_delta(float *dataArr1,int length,int order,float *dataArr2);

// pre_emphasis; coef 0.97
void util_preEmphasis(float *vArr1,int length,float coef,float *vArr2);

// wave相关
int util_readWave(char *name,float **dataArr);
void util_writeWave(char *name,float *dataArr,int length);


#ifdef __cplusplus
}
#endif

#endif