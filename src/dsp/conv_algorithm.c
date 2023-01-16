// clang 

#include <math.h>
#include <string.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "fft_algorithm.h"
#include "conv_algorithm.h"

struct OpaqueConv{
	FFTObj fftObj;

	int fftLength;

	float *dataArr1;
	float *dataArr2;

	float *vRealArr1; // fft相关缓存
	float *vImageArr1;
	float *vRealArr2;
	float *vImageArr2;

};

static int __calFastMethod(int length1,int length2);

static int _conv_direct(float *vArr1,int length1,float *vArr2,int length2, 
				ConvModeType mode,float *vArr3);

static int _convObj_fft(ConvObj convObj,float *vArr1,int length1,float *vArr2,int length2, 
				ConvModeType mode,float *vArr3);

static void _convObj_dealFFT(ConvObj convObj,int fftLength);

int convObj_new(ConvObj *convObj){
	int status=0;
	ConvObj conv=NULL;

	conv=*convObj=(ConvObj )calloc(1, sizeof(struct OpaqueConv ));

	return status;
}

int convObj_conv(ConvObj convObj,float *vArr1,int length1,float *vArr2,int length2,
				ConvModeType *mode,ConvMethodType *method,
				float *vArr3){
	int len=0;
	int fftLength=0;

	ConvMethodType _method=ConvMethod_Auto;
	ConvModeType _mode=ConvMode_Full;

	if(method){
		_method=*method;
	}

	if(mode){
		_mode=*mode;
	}

	if(_method==ConvMethod_Auto){
		fftLength=__calFastMethod(length1,length2);
		if(fftLength){
			_method=ConvMethod_FFT;
		}
		else{
			_method=ConvMethod_Direct;	
		}
	}

	if(_method==ConvMethod_FFT){
		fftLength=util_ceilPowerTwo(2*(length1>length2?length1:length2));
		_convObj_dealFFT(convObj,fftLength);
		len=_convObj_fft(convObj,vArr1,length1,vArr2,length2,_mode,vArr3);
		// printf("fft\n");
	}
	else{ // direct
		len=_conv_direct(vArr1,length1,vArr2,length2,_mode,vArr3);
		// printf("direct\n");
	}

	return len;
}

static void _convObj_dealFFT(ConvObj convObj,int fftLength){
	FFTObj fftObj=NULL;
	int radix2Exp=0;

	float *vRealArr1=NULL; // fft相关缓存
	float *vImageArr1=NULL;
	float *vRealArr2=NULL;
	float *vImageArr2=NULL;

	float *dataArr1=NULL;
	float *dataArr2=NULL;

	if(convObj->fftLength!=fftLength){
		fftObj=convObj->fftObj;

		vRealArr1=convObj->vRealArr1;
		vImageArr1=convObj->vImageArr1;
		vRealArr2=convObj->vRealArr2;
		vImageArr2=convObj->vImageArr2;

		dataArr1=convObj->dataArr1;
		dataArr2=convObj->dataArr2;

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
		
		convObj->fftObj=fftObj;
		convObj->fftLength=fftLength;
		// printf("fftLength is %d\n",fftLength);

		convObj->vRealArr1=vRealArr1;
		convObj->vImageArr1=vImageArr1;
		convObj->vRealArr2=vRealArr2;
		convObj->vImageArr2=vImageArr2;

		convObj->dataArr1=dataArr1;
		convObj->dataArr2=dataArr2;
	}
}

// 卷积定理
static int _convObj_fft(ConvObj convObj,float *vArr1,int length1,float *vArr2,int length2, 
				ConvModeType mode,float *vArr3){
	FFTObj fftObj=NULL;
	int fftLength=0;

	int len=0;

	float *vRealArr1=NULL; // fft相关缓存
	float *vImageArr1=NULL;
	float *vRealArr2=NULL;
	float *vImageArr2=NULL;

	float *dataArr1=NULL;
	float *dataArr2=NULL;

	int start=0;

	fftObj=convObj->fftObj;
	fftLength=convObj->fftLength;

	vRealArr1=convObj->vRealArr1;
	vImageArr1=convObj->vImageArr1;
	vRealArr2=convObj->vRealArr2;
	vImageArr2=convObj->vImageArr2;

	dataArr1=convObj->dataArr1;
	dataArr2=convObj->dataArr2;

	// ifft(fft(A)*fft(B))
	memcpy(dataArr1, vArr1, sizeof(float )*length1);
	memcpy(dataArr2, vArr2, sizeof(float )*length2);
	fftObj_fft(fftObj, dataArr1, NULL, vRealArr1, vImageArr1);
	fftObj_fft(fftObj, dataArr2, NULL, vRealArr2, vImageArr2);

	__vcmul(vRealArr1, vImageArr1, vRealArr2, vImageArr2, fftLength, vRealArr2, vImageArr2);

	fftObj_ifft(fftObj, vRealArr2, vImageArr2, vRealArr1, vImageArr1);

	if(mode==ConvMode_Full){
		len=length1+length2-1;
		memcpy(vArr3, vRealArr1, sizeof(float )*len);
	}
	else{
		len=length1;
		start=length2/2-(length2&1?0:1);
		if(mode==ConvMode_Valid){
			len=length1-length2+1;
			len=(len>0?len:0);
			start=0;
		}

		if(len){
			memcpy(vArr3, vRealArr1+(length2-1-start), sizeof(float )*len);
		}
	}

	return len;
}

/***
	关于mode 以vArr1为准 N/M
	full N+M-1
	same N
	valid N-M+1
****/
static int _conv_direct(float *vArr1,int length1,float *vArr2,int length2, 
				ConvModeType mode,float *vArr3){
	int len=0;
	int start=0;

	if(mode==ConvMode_Full){
		len=length1+length2-1;
		start=-length2+1;
		for(int i=start,n=0;i<length1;i++,n++){
			for(int j=length2-1,k=i;j>=0;j--,k++){
				float _value=0;

				if(k>=0&&k<length1){
					_value=vArr1[k];
				}
				vArr3[n]+=_value*vArr2[j];
			}
		}
	}
	else{
		len=length1;
		start=length2/2-(length2&1?0:1);
		if(mode==ConvMode_Valid){
			len=length1-length2+1;
			start=0;
		}

		
		for(int i=-start,n=0;i<len-start;i++,n++){
			for(int j=length2-1,k=i;j>=0;j--,k++){
				float _value=0;

				if(k>=0&&k<length1){
					_value=vArr1[k];
				}
				vArr3[n]+=_value*vArr2[j];
			}
		}
	}

	return len;
}

/***
	direct N*M
	fft N*logN*3
****/
static int __calFastMethod(int length1,int length2){
	int fftLength=0;

	long long dirTime=0;
	long long fftTime=0;

	fftLength=(length1>length2?length1:length2);
	fftLength=util_ceilPowerTwo(2*fftLength);

	dirTime=length1*length2;
	fftTime=3*fftLength*log2f(fftLength)+fftLength;

	if(fftTime>=dirTime){
		fftLength=0;
	}

	return fftLength;
}

void convObj_free(ConvObj convObj){
	FFTObj fftObj=NULL;

	float *vRealArr1=NULL; // fft相关缓存
	float *vImageArr1=NULL;
	float *vRealArr2=NULL;
	float *vImageArr2=NULL;

	float *dataArr1=NULL;
	float *dataArr2=NULL;

	if(!convObj){
		return;
	}

	fftObj=convObj->fftObj;

	vRealArr1=convObj->vRealArr1;
	vImageArr1=convObj->vImageArr1;
	vRealArr2=convObj->vRealArr2;
	vImageArr2=convObj->vImageArr2;

	dataArr1=convObj->dataArr1;
	dataArr2=convObj->dataArr2;

	fftObj_free(fftObj);

	free(vRealArr1);
	free(vImageArr1);
	free(vRealArr2);
	free(vImageArr2);

	free(dataArr1);
	free(dataArr2);

	free(convObj);
}

