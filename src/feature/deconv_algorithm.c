// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_correct.h"
#include "../dsp/fft_algorithm.h"

#include "deconv_algorithm.h"

struct OpaqueDeconv{
	int num;
	int timeLength;

	FFTObj devFFTObj;
	int devFFTLength;

	float *devDataArr; // spectral mag&fft mag

	float *devRealArr1; // fft
	float *devImageArr1;

	float *devRealArr2; // ifft
	float *devImageArr2;

};

static void __deconvObj_dealDeconv(DeconvObj deconvObj);

int deconvObj_new(DeconvObj *deconvObj,int num){
	int status=0;

	FFTObj devFFTObj=NULL;
	int devFFTLength=0;

	float *devDataArr=NULL; // spectral mag&fft mag

	float *devRealArr1=NULL; // fft
	float *devImageArr1=NULL;

	float *devRealArr2=NULL; // ifft
	float *devImageArr2=NULL;

	int radix2Exp=0;

	if(num<2){
		printf("num is error!!!\n");
		return -1;
	}

	DeconvObj dec=NULL;

	dec=*deconvObj=(DeconvObj )calloc(1, sizeof(struct OpaqueDeconv ));

	devFFTLength=util_ceilPowerTwo(2*num);

	radix2Exp=util_powerTwoBit(devFFTLength);
	fftObj_new(&devFFTObj, radix2Exp);

	devDataArr=__vnew(devFFTLength, NULL);

	devRealArr1=__vnew(devFFTLength, NULL);
	devImageArr1=__vnew(devFFTLength, NULL);

	devRealArr2=__vnew(devFFTLength, NULL);
	devImageArr2=__vnew(devFFTLength, NULL);

	dec->num=num;

	dec->devFFTObj=devFFTObj;
	dec->devFFTLength=devFFTLength;

	dec->devDataArr=devDataArr;

	dec->devRealArr1=devRealArr1;
	dec->devImageArr1=devImageArr1;

	dec->devRealArr2=devRealArr2;
	dec->devImageArr2=devImageArr2;

	return status;
}

void deconvObj_setTimeLength(DeconvObj deconvObj,int timeLength){

	deconvObj->timeLength=timeLength;
}

/***
	mDataArr1 mag/power
	mDataArr2 timbre(formant) 
	mDataArr3 pitch
****/
void deconvObj_deconv(DeconvObj deconvObj,float *mDataArr1,float *mDataArr2,float *mDataArr3){
	int timeLength=0;
	int num=0;
	
	FFTObj devFFTObj=NULL;
	int devFFTLength=0;

	float *devDataArr=NULL; // cqt mag&fft mag

	float *devRealArr1=NULL; // fft
	float *devImageArr1=NULL;

	float *devRealArr2=NULL; // ifft
	float *devImageArr2=NULL;

	timeLength=deconvObj->timeLength;
	num=deconvObj->num;

	// 1. deal cache
	__deconvObj_dealDeconv(deconvObj);
	devFFTObj=deconvObj->devFFTObj;
	devFFTLength=deconvObj->devFFTLength;

	devDataArr=deconvObj->devDataArr;

	devRealArr1=deconvObj->devRealArr1;
	devImageArr1=deconvObj->devImageArr1;

	devRealArr2=deconvObj->devRealArr2;
	devImageArr2=deconvObj->devImageArr2;

	for(int i=0;i<timeLength;i++){
		// 2. fft&mag
		memset(devDataArr, 0, sizeof(float )*devFFTLength);
		memcpy(devDataArr, mDataArr1+i*num, sizeof(float )*num);

		fftObj_fft(devFFTObj, devDataArr, NULL, devRealArr1, devImageArr1);
		__vcabs(devRealArr1, devImageArr1, devFFTLength, devDataArr);

		// 3. timbre --> real(ifft)
		fftObj_ifft(devFFTObj, devDataArr, NULL, devRealArr2, devImageArr2);
		memcpy(mDataArr2+i*num, devRealArr2, sizeof(float )*num);

		// 4. pitch --> real(ifft)
		for(int j=0;j<devFFTLength;j++){
			float _value=0;

			_value=devDataArr[j];
			if(_value<1e-16){
				_value=1e-16;
			}

			// devRealArr1[j]/=(devDataArr[j]+1e-16);
			// devImageArr1[j]/=(devDataArr[j]+1e-16);
			devRealArr1[j]/=_value;
			devImageArr1[j]/=_value;
		}

		fftObj_ifft(devFFTObj, devRealArr1, devImageArr1, devRealArr2, devImageArr2);
		memcpy(mDataArr3+i*num, devRealArr2, sizeof(float )*num);
	}
}

static void __deconvObj_dealDeconv(DeconvObj deconvObj){
	FFTObj devFFTObj=NULL;
	int devFFTLength=0;

	float *devDataArr=NULL;

	float *devRealArr1=NULL;
	float *devImageArr1=NULL;

	float *devRealArr2=NULL;
	float *devImageArr2=NULL;

	int num=0;
	int radix2Exp=0;

	num=deconvObj->num;
	devFFTLength=util_ceilPowerTwo(2*num);
	if(devFFTLength!=deconvObj->devFFTLength){
		devFFTObj=deconvObj->devFFTObj;

		devDataArr=deconvObj->devDataArr;

		devRealArr1=deconvObj->devRealArr1;
		devImageArr1=deconvObj->devImageArr1;

		devRealArr2=deconvObj->devRealArr2;
		devImageArr2=deconvObj->devImageArr2;

		fftObj_free(devFFTObj);

		free(devDataArr);

		free(devRealArr1);
		free(devImageArr1);

		free(devRealArr2);
		free(devImageArr2);

		radix2Exp=util_powerTwoBit(devFFTLength);
		fftObj_new(&devFFTObj, radix2Exp);

		devDataArr=__vnew(devFFTLength, NULL);

		devRealArr1=__vnew(devFFTLength, NULL);
		devImageArr1=__vnew(devFFTLength, NULL);

		devRealArr2=__vnew(devFFTLength, NULL);
		devImageArr2=__vnew(devFFTLength, NULL);

		deconvObj->devFFTObj=devFFTObj;
		deconvObj->devFFTLength=devFFTLength;

		deconvObj->devDataArr=devDataArr;

		deconvObj->devRealArr1=devRealArr1;
		deconvObj->devImageArr1=devImageArr1;

		deconvObj->devRealArr2=devRealArr2;
		deconvObj->devImageArr2=devImageArr2;
	}
}

void deconvObj_free(DeconvObj deconvObj){
	FFTObj devFFTObj=NULL;
	float *devDataArr=NULL; // spectral mag&fft mag

	float *devRealArr1=NULL; // fft
	float *devImageArr1=NULL;

	float *devRealArr2=NULL; // ifft
	float *devImageArr2=NULL;

	if(deconvObj){
		devFFTObj=deconvObj->devFFTObj;
		devDataArr=deconvObj->devDataArr;

		devRealArr1=deconvObj->devRealArr1;
		devImageArr1=deconvObj->devImageArr1;

		devRealArr2=deconvObj->devRealArr2;
		devImageArr2=deconvObj->devImageArr2;

		fftObj_free(devFFTObj);
		free(devDataArr);

		free(devRealArr1);
		free(devImageArr1);

		free(devRealArr2);
		free(devImageArr2);

		free(deconvObj);
	}
}












