// clang 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/flux_window.h"
#include "dsp/fft_algorithm.h"

#include "st_algorithm.h"

struct OpaqueST{
	FFTObj fftObj;

	int fftLength;

	int *binArr;
	int binLength;

	float *realArr1; // 2*fftLength
	float *imageArr1;

	float *realArr2; // 缓存ifft fftLength
	float *imageArr2;

	float *mWinArr; // gauss freNum*winLength -> (fftLength/2+1)*fftLength

	// guass params
	float factor; // λ=1 p=1
	float norm;

};

static void _stObj_initWinData(STObj stObj,int fftLength,float factor,float norm);

int stObj_new(STObj *stObj,int radix2Exp,int minIndex,int maxIndex,float *factor,float *norm){
	int status=0;

	STObj st=NULL;

	FFTObj fftObj=NULL;

	int fftLength=0;

	int *binArr=NULL;
	int binLength=0;

	float *realArr1=NULL; // 2*fftLength
	float *imageArr1=NULL;

	float *realArr2=NULL; // 缓存ifft fftLength
	float *imageArr2=NULL;

	float _factor=1; // λ=1 p=1
	float _norm=1;

	st=*stObj=(STObj )calloc(1, sizeof(struct OpaqueST ));

	if(factor){
		if(*factor>0){
			_factor=*factor;
		}
	}

	if(norm){
		if(*norm>0){
			_norm=*norm;
		}
	}
	
	fftLength=(1<<radix2Exp);
	fftObj_new(&fftObj, radix2Exp);

	realArr1=__vnew(2*fftLength, NULL);
	imageArr1=__vnew(2*fftLength, NULL);

	realArr2=__vnew(fftLength, NULL);
	imageArr2=__vnew(fftLength, NULL);

	binArr=__vnewi(fftLength, NULL);

	_stObj_initWinData(st,fftLength,_factor,_norm);

	// deal binArr
	if(minIndex>=maxIndex||minIndex<0||maxIndex>fftLength/2){
		minIndex=0;
		maxIndex=fftLength/2;
	}

	for(int i=minIndex,j=0;i<=maxIndex;i++,j++){
		binArr[j]=i;
	}
	binLength=maxIndex-minIndex+1;

	st->binArr=binArr;
	st->binLength=binLength;

	st->fftObj=fftObj;
	st->fftLength=fftLength;

	st->realArr1=realArr1;
	st->imageArr1=imageArr1;

	st->realArr2=realArr2;
	st->imageArr2=imageArr2;

	return status;
}

void stObj_useBinArr(STObj stObj,int *binArr,int length){
	int fftLength=0;
	int flag=1;

	fftLength=stObj->fftLength;
	for(int i=0;i<length;i++){
		if(binArr[i]>fftLength/2||binArr[i]<0){
			flag=0;
			break;
		}
	}

	if(flag){
		memcpy(stObj->binArr, binArr, sizeof(int )*length);
		stObj->binLength=length;
	}
}

// update mWinArr
void stObj_setValue(STObj stObj,float factor,float norm){
	float *mWinArr=NULL; 

	if(stObj->factor==factor&&stObj->norm==norm){
		return;
	}

	mWinArr=stObj->mWinArr;
	free(mWinArr);

	_stObj_initWinData(stObj,stObj->fftLength,factor,norm);

}

/***
	1. fft -> outData
	2. loop ifft(outData*winArr)
****/
void stObj_st(STObj stObj,float *dataArr,float *mRealArr,float *mImageArr){
	FFTObj fftObj=NULL;
	int fftLength=0;

	int *binArr=NULL;
	int binLength=0;

	float *realArr1=NULL; // 2*fftLength
	float *imageArr1=NULL;

	float *realArr2=NULL; // fftLength
	float *imageArr2=NULL;

	float *mWinArr=NULL;

	fftObj=stObj->fftObj;
	fftLength=stObj->fftLength;

	realArr1=stObj->realArr1;
	imageArr1=stObj->imageArr1;

	realArr2=stObj->realArr2;
	imageArr2=stObj->imageArr2;

	binArr=stObj->binArr;
	binLength=stObj->binLength;

	mWinArr=stObj->mWinArr;
	
	// 1. fft
	// memset(realArr1, 0, sizeof(float )*2*fftLength);
	// memset(imageArr1, 0, sizeof(float )*2*fftLength);
	fftObj_fft(fftObj, dataArr, NULL, realArr1, imageArr1);
	memcpy(realArr1+fftLength, realArr1, sizeof(float )*fftLength);
	memcpy(imageArr1+fftLength, imageArr1, sizeof(float )*fftLength);

	// 2. loop ifft
	for(int i=0;i<binLength;i++){
		int _index=0;

		_index=binArr[i];
		if(_index!=0){
			memset(imageArr2, 0, sizeof(float )*fftLength); // image 0

			__vcmul(realArr1+_index, imageArr1+_index, mWinArr+_index*fftLength, imageArr2, fftLength, realArr2, imageArr2);
			fftObj_ifft(fftObj, realArr2, imageArr2, mRealArr+i*fftLength, mImageArr+i*fftLength);
		}
		else{
			float _mean=0;

			_mean=__vmean(dataArr, fftLength);
			for(int j=0;j<fftLength;j++){
				mRealArr[i*fftLength+j]=_mean;
			}
		}
	}
}

// gauss data fft之后的数据
static void _stObj_initWinData(STObj stObj,int fftLength,float factor,float norm){
	float *mWinArr=NULL;

	float *eArr1=NULL;
	float *eArr2=NULL;

	float *arr1=NULL;
	float *arr2=NULL;

	float value=0;

	eArr1=__vnew(fftLength, NULL);
	eArr2=__vnew(fftLength, NULL);

	__varange(0, fftLength, 1, &arr1); 
	__varange(-fftLength, 0, 1, &arr2);
	for(int i=0;i<fftLength;i++){
		arr1[i]*=arr1[i];
		arr2[i]*=arr2[i];
	}
	
	// mWinArr=__vnew((fftLength/2+1)*fftLength, NULL);
	mWinArr=(float *)calloc((long long )(fftLength/2+1)*fftLength, sizeof(float ));
	for(int i=1;i<=fftLength/2;i++){ // i=0 不计算
		value=-factor*2*M_PI*M_PI/(powf(i, 2*norm));

		__vmul_value(arr1, value, fftLength, eArr1);
		__vmul_value(arr2, value, fftLength, eArr2);

		__vexp(eArr1, fftLength, eArr1);
		__vexp(eArr2, fftLength, eArr2);

		__vadd(eArr1, eArr2, fftLength, mWinArr+(long long)i*fftLength);
	}

	free(arr1);
	free(arr2);

	free(eArr1);
	free(eArr2);

	stObj->factor=factor;
	stObj->norm=norm;

	stObj->mWinArr=mWinArr;
}

void stObj_free(STObj stObj){
	FFTObj fftObj=NULL;
	int *binArr=NULL;

	float *realArr1=NULL; // 2*fftLength
	float *imageArr1=NULL;

	float *realArr2=NULL; // 缓存ifft fftLength
	float *imageArr2=NULL;

	float *mWinArr=NULL; 

	if(!stObj){
		return;
	}

	fftObj=stObj->fftObj;
	binArr=stObj->binArr;

	realArr1=stObj->realArr1;
	imageArr1=stObj->imageArr1;

	realArr2=stObj->realArr2;
	imageArr2=stObj->imageArr2;

	mWinArr=stObj->mWinArr;

	fftObj_free(fftObj);
	free(binArr);

	free(realArr1);
	free(imageArr1);

	free(realArr2);
	free(imageArr2);

	free(mWinArr);

	free(stObj);
}











