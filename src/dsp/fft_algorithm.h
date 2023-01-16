// clang -g 

#ifndef FFT_ALGORITHM_H
#define FFT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef enum{
	FFTParamError=-100,
	FFTExecError=-101,

} FFTError;

typedef struct OpaqueFFT *FFTObj;

int fftObj_new(FFTObj *fftObj,int radix2Exp);

int fftObj_getFFTLength(FFTObj fftObj);

void fftObj_fft(FFTObj fftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2);
void fftObj_ifft(FFTObj fftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2);
void fftObj_dct(FFTObj fftObj,float *dataArr1,float *dataArr2,int isNorm);
void fftObj_idct(FFTObj fftObj,float *dataArr1,float *dataArr2,int isNorm);

void fftObj_free(FFTObj fftObj);
void fftObj_debug(FFTObj fftObj);

#ifdef __cplusplus
}
#endif

#endif
