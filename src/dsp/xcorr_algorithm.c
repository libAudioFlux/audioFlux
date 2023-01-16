// clang 

#include <math.h>
#include <string.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "fft_algorithm.h"
#include "xcorr_algorithm.h"

struct OpaqueXcorr{
	FFTObj fftObj;

	int fftLength;

	float *dataArr1;
	float *dataArr2;

	float *vRealArr1; // fft相关缓存
	float *vImageArr1;
	float *vRealArr2;
	float *vImageArr2;

};

static int __calFastMethod(int length);

static void _xcorrObj_direct(XcorrObj xcorrObj,float *vArr1,float *vArr2,int length, 
						float *vArr3);

static void _xcorrObj_fft(XcorrObj xcorrObj,float *vArr1,float *vArr2,int length, 
						float *vArr3);

static void _xcorrObj_dealFFT(XcorrObj xcorrObj,int fftLength);

int xcorrObj_new(XcorrObj *xcorrObj){
	int status=0;
	XcorrObj xcorr=NULL;

	xcorr=*xcorrObj=(XcorrObj )calloc(1, sizeof(struct OpaqueXcorr ));

	return status;
}

int xcorrObj_xcorr(XcorrObj xcorrbj,float *vArr1,float *vArr2,int length,
				XcorrNormalType *normType,
				float *vArr3,float *maxValue){
	int index=0;
	
	int fftLength=0;
	XcorrNormalType _normType=XcorrNormal_Coeff;

	if(normType){
		_normType=*normType;
	}

	// fftLength=__calFastMethod(length);
	fftLength=util_ceilPowerTwo(2*length);
	_xcorrObj_dealFFT(xcorrbj,fftLength);

	if(vArr1){ 
		_xcorrObj_fft(xcorrbj,vArr1,vArr2,length,vArr3);
		if(_normType==XcorrNormal_Coeff){
			float sum1=0;
			float sum2=0;
			float scale=0;

			float *_vArr1=NULL;
			float *_vArr2=NULL;

			_vArr1=__vnew(length, NULL);
			if(vArr2){
				_vArr2=__vnew(length, NULL);
			}
			
			memcpy(_vArr1, vArr1, sizeof(float )*length);
			if(vArr2){
				memcpy(_vArr2, vArr2, sizeof(float )*length);
			}

			for(int i=0;i<length;i++){
				_vArr1[i]*=_vArr1[i];
			}

			if(vArr2){
				for(int i=0;i<length;i++){
					_vArr2[i]*=_vArr2[i];
				}
			}

			sum1=__vsum(_vArr1, length);
			if(vArr2){
				sum2=__vsum(_vArr2, length);
			}
			else{
				sum2=sum1;
			}

			scale=sqrtf(sum1*sum2);

			__vdiv_value(vArr3, scale, 2*length-1, NULL);

			free(_vArr1);
			free(_vArr2);
		}

		index=__vmax(vArr3, 2*length-1, maxValue);
	}

	return index;
}

static void _xcorrObj_dealFFT(XcorrObj xcorrbj,int fftLength){
	FFTObj fftObj=NULL;
	int radix2Exp=0;

	float *vRealArr1=NULL; // fft相关缓存
	float *vImageArr1=NULL;
	float *vRealArr2=NULL;
	float *vImageArr2=NULL;

	float *dataArr1=NULL;
	float *dataArr2=NULL;

	if(xcorrbj->fftLength!=fftLength){
		fftObj=xcorrbj->fftObj;

		vRealArr1=xcorrbj->vRealArr1;
		vImageArr1=xcorrbj->vImageArr1;
		vRealArr2=xcorrbj->vRealArr2;
		vImageArr2=xcorrbj->vImageArr2;

		dataArr1=xcorrbj->dataArr1;
		dataArr2=xcorrbj->dataArr2;

		fftObj_free(fftObj);

		free(vRealArr1);
		free(vImageArr1);
		free(vRealArr2);
		free(vImageArr2);

		free(dataArr1);
		free(dataArr2);

		radix2Exp=util_powerTwoBit(fftLength);
		fftObj_new(&fftObj, radix2Exp);

		vRealArr1=__vnew(fftLength, NULL);
		vImageArr1=__vnew(fftLength, NULL);
		vRealArr2=__vnew(fftLength, NULL);
		vImageArr2=__vnew(fftLength, NULL);

		dataArr1=__vnew(fftLength, NULL);
		dataArr2=__vnew(fftLength, NULL);
		
		xcorrbj->fftObj=fftObj;
		xcorrbj->fftLength=fftLength;
		printf("fftLength is %d\n",fftLength);

		xcorrbj->vRealArr1=vRealArr1;
		xcorrbj->vImageArr1=vImageArr1;
		xcorrbj->vRealArr2=vRealArr2;
		xcorrbj->vImageArr2=vImageArr2;

		xcorrbj->dataArr1=dataArr1;
		xcorrbj->dataArr2=dataArr2;
	}
}

// direct 计算
static void _xcorrObj_direct(XcorrObj xcorrObj,float *vArr1,float *vArr2,int length, 
						float *vArr3){

}

// ifft(fft(xn)*conj(fft(hn)))
static void _xcorrObj_fft(XcorrObj xcorrObj,float *vArr1,float *vArr2,int length, 
						float *vArr3){
	FFTObj fftObj=NULL;
	int fftLength=0;

	float *vRealArr1=NULL; // fft相关缓存
	float *vImageArr1=NULL;
	float *vRealArr2=NULL;
	float *vImageArr2=NULL;

	float *dataArr1=NULL;
	float *dataArr2=NULL;

	fftObj=xcorrObj->fftObj;
	fftLength=xcorrObj->fftLength;

	vRealArr1=xcorrObj->vRealArr1;
	vImageArr1=xcorrObj->vImageArr1;
	vRealArr2=xcorrObj->vRealArr2;
	vImageArr2=xcorrObj->vImageArr2;

	dataArr1=xcorrObj->dataArr1;
	dataArr2=xcorrObj->dataArr2;

	// ifft(fft(A)*fft(B))
	memcpy(dataArr1, vArr1, sizeof(float )*length);
	if(vArr2){
		memcpy(dataArr2, vArr2, sizeof(float )*length);
	}
	
	fftObj_fft(fftObj, dataArr1, NULL, vRealArr1, vImageArr1);
	if(vArr2){
		fftObj_fft(fftObj, dataArr2, NULL, vRealArr2, vImageArr2);
	}

	if(vArr2){ // xcorr conj dot
		for(int i=0;i<fftLength;i++){
			vImageArr2[i]=-vImageArr2[i];
		}

		__vcmul(vRealArr1, vImageArr1, vRealArr2, vImageArr2, fftLength, vRealArr2, vImageArr2);
	}
	else{ // autocorr abs^2
		memset(vImageArr2,0,sizeof(float )*fftLength);

		__vcabs(vRealArr1, vImageArr1, fftLength, vRealArr2);
		for(int i=0;i<fftLength;i++){
			vRealArr2[i]*=vRealArr2[i];
		}
	}

	fftObj_ifft(fftObj, vRealArr2, vImageArr2, vRealArr1, vImageArr1);

	// fftLength-(length-1)~fftLength is negative 0~length is positive
	for(int i=fftLength-(length-1),j=0;i<fftLength;i++,j++){
		vArr3[j]=vRealArr1[i];
	}
	for(int i=0;i<length;i++){
		vArr3[length-1+i]=vRealArr1[i];
	}

}

void xcorrObj_free(XcorrObj xcorrObj){
	FFTObj fftObj=NULL;

	float *vRealArr1=NULL; // fft相关缓存
	float *vImageArr1=NULL;
	float *vRealArr2=NULL;
	float *vImageArr2=NULL;

	float *dataArr1=NULL;
	float *dataArr2=NULL;

	if(!xcorrObj){
		return;
	}

	fftObj=xcorrObj->fftObj;

	vRealArr1=xcorrObj->vRealArr1;
	vImageArr1=xcorrObj->vImageArr1;
	vRealArr2=xcorrObj->vRealArr2;
	vImageArr2=xcorrObj->vImageArr2;

	dataArr1=xcorrObj->dataArr1;
	dataArr2=xcorrObj->dataArr2;

	fftObj_free(fftObj);

	free(vRealArr1);
	free(vImageArr1);
	free(vRealArr2);
	free(vImageArr2);

	free(dataArr1);
	free(dataArr2);

	free(xcorrObj);
}

/***
	direct N*N
	fft N*logN*3
****/
static int __calFastMethod(int length){
	int fftLength=0;

	long long dirTime=0;
	long long fftTime=0;

	fftLength=util_ceilPowerTwo(2*length);

	dirTime=length*length;
	fftTime=3*fftLength*log2f(fftLength)+fftLength;

	if(fftTime>=dirTime){
		fftLength=0;
	}

	return fftLength;
}







