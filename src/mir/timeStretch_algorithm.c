// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "../util/flux_util.h"
#include "../dsp/phase_vocoder.h"
#include "../stft_algorithm.h"

#include "timeStretch_algorithm.h"

struct OpaqueTimeStretch{
	STFTObj stftObj;

	int radix2Exp;
	int fftLength;
	int slideLength;
	WindowType windowType;

	int timeLength1;
	int timeLength2;

	float *mRealArr1; // timeLength1
	float *mImageArr1;

	float *mRealArr2; // timeLength2
	float *mImageArr2;

};

int timeStretchObj_new(TimeStretchObj *timeStretchObj,int *radix2Exp,int *slideLength,WindowType *windowType){
	int status=0;

	int fftLength=0;
	int _radix2Exp=12;
	int _slideLength=0;
	WindowType _windowType=Window_Hann;

	STFTObj stftObj=NULL;
	TimeStretchObj ts=NULL;

	ts=*timeStretchObj=(TimeStretchObj )calloc(1,sizeof(struct OpaqueTimeStretch ));

	if(radix2Exp){
		if(*radix2Exp>=1&&*radix2Exp<=30){
			_radix2Exp=*radix2Exp;
		}
	}

	fftLength=1<<_radix2Exp;
	_slideLength=fftLength/4;
	if(slideLength){
		if(*slideLength>0){ // &&*slideLength<=fftLength support not overlap
			_slideLength=*slideLength;
		}
	}

	if(windowType){
		_windowType=*windowType;
	}

	stftObj_new(&stftObj, _radix2Exp, &_windowType, &_slideLength, NULL);

	ts->stftObj=stftObj;

	ts->radix2Exp=_radix2Exp;
	ts->fftLength=fftLength;
	ts->slideLength=_slideLength;
	ts->windowType=_windowType;

	return status;
}

int timeStretchObj_calDataCapacity(TimeStretchObj timeStretchObj,float rate,int dataLength){

	return ceilf(dataLength/rate)+timeStretchObj->fftLength;
}

int timeStretchObj_timeStretch(TimeStretchObj timeStretchObj,float rate,float *dataArr1,int dataLength,float *dataArr2){
	int timeLength1=0;
	int timeLength2=0;

	int fftLength=0;
	int slideLength=0;

	float *mRealArr1=NULL; // timeLength1
	float *mImageArr1=NULL;

	float *mRealArr2=NULL; // timeLength2
	float *mImageArr2=NULL;

	if(rate<0){
		return 0;
	}

	fftLength=timeStretchObj->fftLength;
	slideLength=timeStretchObj->slideLength;

	timeLength1=stftObj_calTimeLength(timeStretchObj->stftObj, dataLength);
	timeLength2=ceilf(timeLength1/rate);

	mRealArr1=timeStretchObj->mRealArr1;
	mImageArr1=timeStretchObj->mImageArr1;

	mRealArr2=timeStretchObj->mRealArr2;
	mImageArr2=timeStretchObj->mImageArr2;

	if(timeStretchObj->timeLength1<timeLength1||
		timeStretchObj->timeLength1>timeLength1*2){
		free(mRealArr1);
		free(mImageArr1);

		mRealArr1=__vnew(timeLength1*fftLength, NULL);
		mImageArr1=__vnew(timeLength1*fftLength, NULL);

		timeStretchObj->timeLength1=timeLength1;

		timeStretchObj->mRealArr1=mRealArr1;
		timeStretchObj->mImageArr1=mImageArr1;
	}

	if(timeStretchObj->timeLength2<timeLength2||
		timeStretchObj->timeLength2>timeLength2*2){
		free(mRealArr2);
		free(mImageArr2);

		mRealArr2=__vnew(timeLength2*fftLength, NULL);
		mImageArr2=__vnew(timeLength2*fftLength, NULL);

		timeStretchObj->timeLength2=timeLength2;

		timeStretchObj->mRealArr2=mRealArr2;
		timeStretchObj->mImageArr2=mImageArr2;
	}
	
	// 1. stft
	stftObj_stft(timeStretchObj->stftObj, dataArr1, dataLength, mRealArr1, mImageArr1);

	// 2. phase vocoder
	phase_vocoder(mRealArr1,mImageArr1,timeLength1,fftLength,slideLength,rate,mRealArr2,mImageArr2);

	// 3. istft
	stftObj_istft(timeStretchObj->stftObj, mRealArr2, mImageArr2, timeLength2, 0, dataArr2);

	return roundf(dataLength/rate);
}

void timeStretchObj_free(TimeStretchObj timeStretchObj){

	if(timeStretchObj){
		stftObj_free(timeStretchObj->stftObj);

		free(timeStretchObj->mRealArr1);
		free(timeStretchObj->mImageArr1);

		free(timeStretchObj->mRealArr2);
		free(timeStretchObj->mImageArr2);

		free(timeStretchObj);
	}
}









