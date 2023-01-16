// clang -g -c 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_complex.h"

#include "flux_window.h"
#include "fft_algorithm.h"
#include "hilbert_algorithm.h"

struct OpaqueHilbert{
	FFTObj fftObj;
	int fftLength; // 2*dataLength

	float *realArr1; // fft result
	float *imageArr1; 

	float *hArr; // h index

};

int hilbertObj_new(HilbertObj *hilbertObj,int radix2Exp){
	int status=0;

	FFTObj fftObj=NULL;
	int fftLength=0;

	float *realArr1=NULL; // fft result
	float *imageArr1=NULL; 

	float *hArr=NULL; // h index

	HilbertObj hilbert=NULL;

	hilbert=*hilbertObj=(HilbertObj )calloc(1, sizeof(struct OpaqueHilbert ));

	status=fftObj_new(&fftObj, radix2Exp);
	fftLength=(1<<radix2Exp);

	realArr1=__vnew(fftLength, NULL);
	imageArr1=__vnew(fftLength, NULL);

	hArr=__vnew(fftLength, NULL);
	hArr[0]=1;
	hArr[fftLength/2]=1;
	for(int i=1;i<fftLength/2;i++){
		hArr[i]=2;
	}

	hilbert->fftObj=fftObj;
	hilbert->fftLength=fftLength;

	hilbert->realArr1=realArr1;
	hilbert->imageArr1=imageArr1;

	hilbert->hArr=hArr;

	return status;
}

void hilbertObj_hilbert(HilbertObj hilbertObj,float *dataArr,
						float *realArr3,float *imageArr3){
	FFTObj fftObj=NULL;
	int fftLength=0; // 2*dataLength

	float *realArr1=NULL; // fft result
	float *imageArr1=NULL; 

	float *hArr=NULL;

	fftObj=hilbertObj->fftObj;
	fftLength=hilbertObj->fftLength;

	realArr1=hilbertObj->realArr1;
	imageArr1=hilbertObj->imageArr1;

	hArr=hilbertObj->hArr;

	// 1. fft
	fftObj_fft(fftObj, dataArr, NULL, realArr1, imageArr1);

	// 2. ifft
	for(int i=0;i<fftLength;i++){
		realArr1[i]*=hArr[i];
		imageArr1[i]*=hArr[i];
	}

	fftObj_ifft(fftObj, realArr1, imageArr1, realArr3, imageArr3);
}

void hilbertObj_free(HilbertObj hilbertObj){
	FFTObj fftObj=NULL;

	float *realArr1=NULL; // fft result
	float *imageArr1=NULL; 

	float *hArr=NULL; 

	if(!hilbertObj){
		return;
	}

	fftObj=hilbertObj->fftObj;

	realArr1=hilbertObj->realArr1;
	imageArr1=hilbertObj->imageArr1;

	hArr=hilbertObj->hArr;

	fftObj_free(fftObj);

	free(realArr1);
	free(imageArr1);

	free(hArr);

	free(hilbertObj);
}









