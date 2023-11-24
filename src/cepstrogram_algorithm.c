// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/fft_algorithm.h"

#include "stft_algorithm.h"
#include "cepstrogram_algorithm.h"

struct OpaqueCepstrogram{
	STFTObj stftObj;
	FFTObj fftObj;

	int fftLength; // fftLength,timeLength
	int timeLength;

	int radix2Exp; // stft
	WindowType windowType;
	int slideLength;

	// stft
	float *mYArr1; // ifft(log(S))
	float *mRealArr1; // stft r,i(timeLength*fftLength)
	float *mImageArr1;
	float *mSArr1; // power 
	
	float *mYArr2; //  envelope
	float *mRealArr2; 
	float *mImageArr2;
	
	float *mYArr3; // details
	float *mRealArr3;
	float *mImageArr3;

	int isDebug;

};

static void __cepstrogramObj_spectrogram(CepstrogramObj cepstrogramObj,int cepNum,float *dataArr,int dataLength,
										float *mRealArr,float *mImageArr,int nLength,
										float *mDataArr1,float *mDataArr2,float *mDataArr3);

/***
	radix2Exp 12 
	samplate 32000
	WindowType "rect"
	slideLength 1024
****/
int cepstrogramObj_new(CepstrogramObj *cepstrogramObj,int radix2Exp,
					WindowType *windowType,int *slideLength){
	int status=0;
	CepstrogramObj ceps=NULL;

	STFTObj stftObj=NULL;
	FFTObj fftObj=NULL;

	int fftLength=0;

	WindowType _windowType=Window_Rect;
	int _slideLength=0;

	if(radix2Exp<1||radix2Exp>30){
		status=-100;
		printf("radix2Exp is error!\n");
		return status;
	}

	fftLength=1<<radix2Exp;

	if(windowType){
		_windowType=*windowType;
	}

	_slideLength=fftLength/4;
	if(slideLength){
		if(*slideLength>0){ // &&*slideLength<=fftLength ->support not overlap
			_slideLength=*slideLength;
		}
	}

	ceps=*cepstrogramObj=(CepstrogramObj )calloc(1, sizeof(struct OpaqueCepstrogram ));

	stftObj_new(&stftObj,radix2Exp,&_windowType,&_slideLength,NULL);
	fftObj_new(&fftObj, radix2Exp);

	ceps->stftObj=stftObj;
	ceps->fftObj=fftObj;

	ceps->fftLength=fftLength;

	ceps->radix2Exp=radix2Exp;
	ceps->windowType=_windowType;
	ceps->slideLength=_slideLength;

	return status;
}

int cepstrogramObj_calTimeLength(CepstrogramObj cepstrogramObj,int dataLength){
	int timeLength=0;

	timeLength=stftObj_calTimeLength(cepstrogramObj->stftObj,dataLength);
	return timeLength;
}

void cepstrogramObj_cepstrogram(CepstrogramObj cepstrogramObj,int cepNum,float *dataArr,int dataLength,
								float *mDataArr1,float *mDataArr2,float *mDataArr3){

	__cepstrogramObj_spectrogram(cepstrogramObj,cepNum,dataArr,dataLength,
								NULL,NULL,0,
								mDataArr1,mDataArr2,mDataArr3);
}

void cepstrogramObj_cepstrogram2(CepstrogramObj cepstrogramObj,int cepNum,float *mRealArr,float *mImageArr,int nLength,
								float *mDataArr1,float *mDataArr2,float *mDataArr3){

	__cepstrogramObj_spectrogram(cepstrogramObj,cepNum,NULL,0,
								mRealArr,mImageArr,nLength,
								mDataArr1,mDataArr2,mDataArr3);
}

static void __cepstrogramObj_spectrogram(CepstrogramObj cepstrogramObj,int cepNum,float *dataArr,int dataLength,
										float *mRealArr,float *mImageArr,int nLength,
										float *mDataArr1,float *mDataArr2,float *mDataArr3){
	STFTObj stftObj=NULL;
	FFTObj fftObj=NULL;

	int fftLength=0; // fftLength,timeLength,num
	int timeLength=0;

	float *mYArr1=NULL; // ifft(log(S))
	float *mRealArr1=NULL; // stft r,i(timeLength*fftLength)
	float *mImageArr1=NULL;
	float *mSArr1=NULL; // power 
	
	float *mYArr2=NULL; //  envelope
	float *mRealArr2=NULL; 
	float *mImageArr2=NULL;
	
	float *mYArr3=NULL; // details
	float *mRealArr3=NULL;
	float *mImageArr3=NULL;

	int sFlag=0;

	stftObj=cepstrogramObj->stftObj;
	fftObj=cepstrogramObj->fftObj;

	fftLength=cepstrogramObj->fftLength;

	if(dataArr&&dataLength){
		sFlag=1;
	}

	if(sFlag){
		timeLength=stftObj_calTimeLength(stftObj,dataLength);
	}
	else{
		timeLength=nLength;
	}

	if(cepstrogramObj->timeLength<timeLength||
		cepstrogramObj->timeLength>timeLength*2){ // 更新缓存
		free(cepstrogramObj->mYArr1);
		free(cepstrogramObj->mRealArr1);
		free(cepstrogramObj->mImageArr1);
		free(cepstrogramObj->mSArr1);

		free(cepstrogramObj->mYArr2);
		free(cepstrogramObj->mRealArr2);
		free(cepstrogramObj->mImageArr2);

		free(cepstrogramObj->mYArr3);
		free(cepstrogramObj->mRealArr3);
		free(cepstrogramObj->mImageArr3);
		
		cepstrogramObj->mYArr1=__vnew(timeLength*fftLength, NULL);
		cepstrogramObj->mRealArr1=__vnew(timeLength*fftLength, NULL);
		cepstrogramObj->mImageArr1=__vnew(timeLength*fftLength, NULL);
		cepstrogramObj->mSArr1=__vnew(timeLength*fftLength, NULL);

		cepstrogramObj->mYArr2=__vnew(timeLength*fftLength, NULL);
		cepstrogramObj->mRealArr2=__vnew(timeLength*fftLength, NULL);
		cepstrogramObj->mImageArr2=__vnew(timeLength*fftLength, NULL);

		cepstrogramObj->mYArr3=__vnew(timeLength*fftLength, NULL);
		cepstrogramObj->mRealArr3=__vnew(timeLength*fftLength, NULL);
		cepstrogramObj->mImageArr3=__vnew(timeLength*fftLength, NULL);
	}

	mYArr1=cepstrogramObj->mYArr1;
	mRealArr1=cepstrogramObj->mRealArr1;
	mImageArr1=cepstrogramObj->mImageArr1;
	mSArr1=cepstrogramObj->mSArr1;

	mYArr2=cepstrogramObj->mYArr2;
	mRealArr2=cepstrogramObj->mRealArr2;
	mImageArr2=cepstrogramObj->mImageArr2;

	mYArr3=cepstrogramObj->mYArr3;
	mRealArr3=cepstrogramObj->mRealArr3;
	mImageArr3=cepstrogramObj->mImageArr3;

	// 1. stft
	if(sFlag){
		stftObj_stft(stftObj,dataArr,dataLength,mRealArr1,mImageArr1);
	}
	else{
		memcpy(mRealArr, mRealArr1, sizeof(float )*nLength*fftLength);
		memcpy(mImageArr, mImageArr1, sizeof(float )*nLength*fftLength);
	}

	// logf
	__vcsquare(mRealArr1,mImageArr1,timeLength*fftLength,mSArr1);
	for(int i=0;i<timeLength*fftLength;i++){
		float r1=0;

		r1=mSArr1[i];
		if(r1<1e-16){
			r1=1e-16;
		}
		
		mSArr1[i]=logf(r1);
	}

	// 2. ifft -> mYArr1, mYArr1 is even
	for(int i=0;i<timeLength;i++){
		fftObj_ifft(fftObj, mSArr1+i*fftLength, NULL, mYArr1+i*fftLength, mImageArr1+i*fftLength);
	}

	if(mDataArr1){
		for(int i=0;i<timeLength;i++){
			memcpy(mDataArr1+i*(fftLength/2+1), mYArr1+i*fftLength, sizeof(float )*(fftLength/2+1));
		}

		if(cepstrogramObj->isDebug){
			printf("mDataArr1 is :\n");
			__mdebug(mDataArr1, timeLength, fftLength/2+1, 1);
			printf("\n");
		}
	}

	// 3. deconv
	if(mDataArr2){ // envelope
		for(int i=0;i<timeLength;i++){
			float *_arr1=NULL;
			float *_arr2=NULL;

			_arr1=mYArr1+i*fftLength;
			_arr2=mYArr2+i*fftLength;

			// [0,cepNum] ... [cepNum-1,1]
			memset(_arr2, 0, sizeof(float )*fftLength);
			memcpy(_arr2,_arr1,sizeof(float)*(cepNum+1));
			for(int j=0;j<cepNum;j++){
				_arr2[fftLength-j-1]=_arr2[j+1];
			}

			fftObj_fft(fftObj, _arr2, NULL, mRealArr2+i*fftLength, mImageArr2+i*fftLength);
			memcpy(mDataArr2+i*(fftLength/2+1), mRealArr2+i*fftLength, sizeof(float )*(fftLength/2+1));
		}

		if(cepstrogramObj->isDebug){
			printf("mDataArr2 is :\n");
			__mdebug(mDataArr2, timeLength, fftLength/2+1, 1);
			printf("\n");
		}
	}

	if(mDataArr3){ // details
		for(int i=0;i<timeLength;i++){
			float *_arr1=NULL;
			float *_arr2=NULL;

			_arr1=mYArr1+i*fftLength;
			_arr2=mYArr3+i*fftLength;

			memset(_arr2, 0, sizeof(float )*fftLength);
			memcpy(_arr2+(cepNum+1),_arr1+(cepNum+1),sizeof(float)*(fftLength-2*cepNum));

			fftObj_fft(fftObj, _arr2, NULL, mRealArr3+i*fftLength, mImageArr3+i*fftLength);
			memcpy(mDataArr3+i*(fftLength/2+1), mRealArr3+i*fftLength, sizeof(float )*(fftLength/2+1));
		}

		if(cepstrogramObj->isDebug){
			printf("mDataArr3 is :\n");
			__mdebug(mDataArr3, timeLength, fftLength/2+1, 1);
			printf("\n");
		}
	}

	cepstrogramObj->timeLength=timeLength;
}

void cepstrogramObj_enableDebug(CepstrogramObj cepstrogramObj,int flag){

	cepstrogramObj->isDebug=flag;
}

void cepstrogramObj_free(CepstrogramObj cepstrogramObj){
	STFTObj stftObj=NULL;
	FFTObj fftObj=NULL;

	// stft
	float *mYArr1=NULL; // ifft(log(S))
	float *mRealArr1=NULL; // stft r,i(timeLength*fftLength)
	float *mImageArr1=NULL;
	float *mSArr1=NULL; // power 
	
	float *mYArr2=NULL; //  envelope
	float *mRealArr2=NULL; 
	float *mImageArr2=NULL;
	
	float *mYArr3=NULL; // details
	float *mRealArr3=NULL;
	float *mImageArr3=NULL;

	if(!cepstrogramObj){
		return;
	}

	stftObj=cepstrogramObj->stftObj;
	fftObj=cepstrogramObj->fftObj;

	mYArr1=cepstrogramObj->mYArr1;
	mRealArr1=cepstrogramObj->mRealArr1;
	mImageArr1=cepstrogramObj->mImageArr1;
	mSArr1=cepstrogramObj->mSArr1;

	mYArr2=cepstrogramObj->mYArr2;
	mRealArr2=cepstrogramObj->mRealArr2;
	mImageArr2=cepstrogramObj->mImageArr2;

	mYArr3=cepstrogramObj->mYArr3;
	mRealArr3=cepstrogramObj->mRealArr3;
	mImageArr3=cepstrogramObj->mImageArr3;

	stftObj_free(stftObj);
	fftObj_free(fftObj);

	free(mYArr1);
	free(mRealArr1);
	free(mImageArr1);
	free(mSArr1);

	free(mYArr2);
	free(mRealArr2);
	free(mImageArr2);

	free(mYArr3);
	free(mRealArr3);
	free(mImageArr3);

	free(cepstrogramObj);
}










