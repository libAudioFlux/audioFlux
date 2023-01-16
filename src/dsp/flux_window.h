

#ifndef FLUX_WINDOW_H
#define FLUX_WINDOW_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../flux_base.h"

// 窗函数相关 flag 0 symmetric对称 1 periodic周期 针对fft length extend 1 sample
float *window_createHann(int length,int flag);
float *window_createHamm(int length,int flag);

float *window_createBlackman(int length,int flag);
// beta 5.0
float *window_createKaiser(int length,int flag,float *beta);

// only symmetric
float *window_createBartlett(int length,int flag);
float *window_createTriang(int length,int flag);

float *window_createFlattop(int length,int flag);
// alpha 2.5
float *window_createGauss(int length,int flag,float *alpha);

// 4-term blackman harris~nuttall
float *window_createBlackmanHarris(int length,int flag);
float *window_createBlackmanNuttall(int length,int flag);
// like Bartlett,hann,hamm; only symmetric
float *window_createBartlettHann(int length,int flag);

// only symmetric
float *window_createBohman(int length,int flag);

// alpha>=0&&<=1 0.5 rect~hann ;only symmetric
float *window_createTukey(int length,int flag,float *alpha);

// order相关 
// w1 ws/wp w2 ws/wp atten dB
void window_calKaiserOrder(float w1,float w2,float atten,int *order,float *beta);

float *window_calFFTWindow(WindowType type,int length);


#ifdef __cplusplus
}
#endif

#endif