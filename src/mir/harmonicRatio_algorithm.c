// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"

#include "harmonicRatio_algorithm.h"

struct OpaqueHarmonicRatio{
	FFTObj fftObj;

	int samplate;
	int fftLength; 

	int windowLength; // fftLength/2
	int slideLength;

	int maxLength; // samplate*lowFre
	int minLength; // first zero crossing

	int timeLength;

	float *winDataArr; // windowLength

	// fftLength 
	float *curDataArr;

	float *vArr1; 
	float *vArr2;

	float *realArr1;
	float *imageArr1;

	float *realArr2;
	float *imageArr2;

};

/***
	samplate 32000
	radix2Exp 12 
	windowType hamm
****/
int harmonicRatioObj_new(HarmonicRatioObj *harmonicRatioObj,
						int *samplate,float *lowFre,
						int *radix2Exp,WindowType *windowType,int *slideLength){
	int status=0;

	int _samplate=32000;
	float _lowFre=25;
	int _radix2Exp=12;
	WindowType _windowType=Window_Hamm;
	int _slideLength=0;

	int fftLength=0;
	int maxLength=0;
	int windowLength=0;

	float *winDataArr=NULL;

	float *vArr1=NULL;
	float *vArr2=NULL;

	float *curDataArr=NULL;

	float *realArr1=NULL;
	float *imageArr1=NULL;

	float *realArr2=NULL;
	float *imageArr2=NULL;

	FFTObj fftObj=NULL;
	HarmonicRatioObj hr=NULL;

	hr=*harmonicRatioObj=(HarmonicRatioObj )calloc(1, sizeof(struct OpaqueHarmonicRatio ));

	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
		}
	}

	if(lowFre){
		if(*lowFre>0&&*lowFre<_samplate/2){
			_lowFre=*lowFre;
		}
	}

	if(radix2Exp){
		if(*radix2Exp+1>=1&&*radix2Exp+1<=30){
			_radix2Exp=*radix2Exp+1;
		}
	}

	fftLength=1<<_radix2Exp;
	windowLength=fftLength/2;
	_slideLength=windowLength/4;
	if(slideLength){
		if(*slideLength>0){ // &&*slideLength<=fftLength support not overlap
			_slideLength=*slideLength;
		}
	}

	maxLength=floorf(_samplate/_lowFre);
	if(maxLength>windowLength-1){
		maxLength=windowLength-1;
	}

	fftObj_new(&fftObj, _radix2Exp);

	winDataArr=window_calFFTWindow(_windowType,windowLength);

	vArr1=__vnew(fftLength, NULL);
	vArr2=__vnew(fftLength, NULL);

	curDataArr=__vnew(fftLength, NULL);

	realArr1=__vnew(fftLength, NULL);
	imageArr1=__vnew(fftLength, NULL);

	realArr2=__vnew(fftLength, NULL);
	imageArr2=__vnew(fftLength, NULL);

	hr->fftObj=fftObj;

	hr->samplate=_samplate;
	hr->maxLength=maxLength;

	hr->fftLength=fftLength;
	hr->windowLength=windowLength;
	hr->slideLength=_slideLength;

	hr->winDataArr=winDataArr;

	hr->vArr1=vArr1;
	hr->vArr2=vArr2;

	hr->curDataArr=curDataArr;

	hr->realArr1=realArr1;
	hr->imageArr1=imageArr1;

	hr->realArr2=realArr2;
	hr->imageArr2=imageArr2;

	return status;
}

int harmonicRatioObj_calTimeLength(HarmonicRatioObj harmonicRatioObj,int dataLength){
	int windowLength=0; 
	int slideLength=0;
	int timeLength=0;

	windowLength=harmonicRatioObj->windowLength;
	slideLength=harmonicRatioObj->slideLength;
	if(dataLength<windowLength){
		return 0;
	}

	timeLength=(dataLength-windowLength)/slideLength+1;
	return timeLength;
}

void harmonicRatioObj_harmonicRatio(HarmonicRatioObj harmonicRatioObj,float *dataArr,int dataLength,float *valueArr){
	FFTObj fftObj;

	int fftLength=0; 

	int windowLength=0;
	int slideLength=0;
	int maxLength=0;

	// fftLength 
	float *winDataArr=NULL;

	float *curDataArr=NULL;

	float *vArr1=NULL;
	float *vArr2=NULL;

	float *realArr1=NULL;
	float *imageArr1=NULL;

	float *realArr2=NULL;
	float *imageArr2=NULL;

	int timeLength=0;
	int minIndex=0;

	fftObj=harmonicRatioObj->fftObj;

	fftLength=harmonicRatioObj->fftLength;

	windowLength=harmonicRatioObj->windowLength;
	slideLength=harmonicRatioObj->slideLength;
	maxLength=harmonicRatioObj->maxLength;

	winDataArr=harmonicRatioObj->winDataArr;

	curDataArr=harmonicRatioObj->curDataArr;

	vArr1=harmonicRatioObj->vArr1;
	vArr2=harmonicRatioObj->vArr2;

	realArr1=harmonicRatioObj->realArr1;
	imageArr1=harmonicRatioObj->imageArr1;

	realArr2=harmonicRatioObj->realArr2;
	imageArr2=harmonicRatioObj->imageArr2;

	if(dataLength<windowLength){
		return ;
	}

	timeLength=(dataLength-windowLength)/slideLength+1;
	for(int i=0;i<timeLength;i++){
		float cur=0;

		int index=0;
		float value1=0;
		float value2=0;
		float value3=0;

		// 0. reset
		memset(realArr1, 0, sizeof(float )*fftLength);
		memset(imageArr1, 0, sizeof(float )*fftLength);

		memset(realArr2, 0, sizeof(float )*fftLength);
		memset(imageArr2, 0, sizeof(float )*fftLength);

		// 1. auto corr --> realArr2
		__vmul(dataArr+i*slideLength, winDataArr, windowLength, curDataArr);

		fftObj_fft(fftObj, curDataArr, NULL, realArr1, imageArr1);
		for(int j=0;j<fftLength;j++){
			vArr1[j]=realArr1[j]*realArr1[j]+imageArr1[j]*imageArr1[j];
		}
		fftObj_ifft(fftObj, vArr1, NULL, realArr2, imageArr2);

		// 2. cumsum power -->vArr2
		for(int j=0;j<windowLength;j++){
			cur+=curDataArr[j]*curDataArr[j];
			imageArr2[j]=cur;
		}

		for(int j=windowLength-2,k=0;j>windowLength-maxLength-2;j--,k++){
			vArr2[k]=imageArr2[j];
		}

		// 3. minIndex
		for(int j=2;j<=maxLength;j++){
			if((realArr2[j]>=0&&realArr2[j-1]<=0)||
				(realArr2[j]<=0&&realArr2[j-1]>=0)){

				minIndex=j-1;
				break;
			}
		}

		// 4. gamma -->vArr1
		// printf("gamma=%d \n",i);
		for(int j=minIndex+1,k=0;j<maxLength;j++,k++){
			vArr1[k]=realArr2[j]/sqrtf(realArr2[0]*vArr2[j]+1e-16);
			// printf("%d, %f\n",j,sqrtf(realArr2[0]*vArr2[j]));
		}
		// printf("\n");

		// 5. parab inter
		index=__vmax(vArr1, maxLength-minIndex-1, &value2);
		if(index==0||index==maxLength-minIndex-2){
			valueArr[i]=value2;
		}
		else{
			value1=vArr1[index-1];
			value3=vArr1[index+1];
			util_qaudInterp(value1,value2,value3,valueArr+i);
		}
	}
}

void harmonicRatioObj_free(HarmonicRatioObj harmonicRatioObj){

	if(harmonicRatioObj){
		fftObj_free(harmonicRatioObj->fftObj);

		free(harmonicRatioObj->winDataArr);
		free(harmonicRatioObj->curDataArr);

		free(harmonicRatioObj->vArr1);
		free(harmonicRatioObj->vArr2);

		free(harmonicRatioObj->realArr1);
		free(harmonicRatioObj->imageArr1);

		free(harmonicRatioObj->realArr2);
		free(harmonicRatioObj->imageArr2);

		free(harmonicRatioObj);
	}
}










