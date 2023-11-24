// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"

#include "_pitch_lhs.h"

struct OpaquePitchLHS{
	int isContinue;

	FFTObj fftObj;

	int fftLength; 
	int slideLength;
	int radix2Exp; // fftLength

	int interpFFTLength; // samplate
	int timeLength;

	int minIndex; // edge
	int maxIndex;

	int harmonicCount;
	float *winDataArr; // fftLength

	float *mDbArr; // timeLength*interpFFTLength
	float *mSumArr; 

	// cache data
	float *realArr1; // interpFFTLength
	float *imageArr1;

	float *dataArr1; 

	// continue
	float *tailDataArr; // fftLength
	int tailDataLength;

	float *curDataArr;
	int curDataLength;

	int samplate;
	WindowType winType;

	float lowFre;
	float highFre;

	int isDebug;
};

static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int *timeLength,int *tailLength);

static void __pitchLHSObj_initData(PitchLHSObj pitchLHSObj);
static int __pitchLHSObj_dealData(PitchLHSObj pitchLHSObj,float *dataArr,int dataLength);

static void __pitchLHSObj_calDb(PitchLHSObj pitchLHSObj);
static void __pitchLHSObj_calSum(PitchLHSObj pitchLHSObj);
static void __pitchLHSObj_dealResult(PitchLHSObj pitchLHSObj,float *freArr);

/***
	samplate 32000
	lowFre 32,
	highFre 2000

	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	
	harmonicCount 5 >0
	isContinue 0
****/
int pitchLHSObj_new(PitchLHSObj *pitchLHSObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				int *harmonicCount,
				int *isContinue){
	int status=0;
	
	float _lowFre=32; // 27.5
	float _highFre=2000; // 2093/4186

	int _samplate=32000;
	int _radix2Exp=12;
	int _slideLength=0;
	WindowType _winType=Window_Hamm;
	int _isContinue=0;

	int fftLength=0;
	int interpFFTLength=0;
	int radix2Exp2=0;

	int _hc=5;

	FFTObj fftObj=NULL;
	PitchLHSObj pitch=NULL;

	pitch=*pitchLHSObj=(PitchLHSObj )calloc(1,sizeof(struct OpaquePitchLHS ));

	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
		}
	}

	if(lowFre){
		if(*lowFre>=27){
			_lowFre=*lowFre;
		}
	}

	if(highFre){
		if(*highFre>_lowFre&&*highFre<_samplate/2){
			_highFre=*highFre;
		}
		else{
			_lowFre=32;
			_highFre=2000;
		}
	}

	if(radix2Exp){
		if(*radix2Exp>=1&&*radix2Exp<=30){
			_radix2Exp=*radix2Exp;
		}
	}

	if(harmonicCount){
		if(*harmonicCount>0){
			_hc=*harmonicCount;
		}
	}

	if(windowType){
		_winType=*windowType;
	}

	fftLength=1<<_radix2Exp;
	_slideLength=fftLength/4;
	if(slideLength){
		if(*slideLength>0){ // &&*slideLength<=fftLength support not overlap
			_slideLength=*slideLength;
		}
	}

	if(isContinue){
		_isContinue=*isContinue;
	}

	interpFFTLength=util_roundPowerTwo(_samplate);
	radix2Exp2=util_powerTwoBit(interpFFTLength);

	fftObj_new(&fftObj, radix2Exp2);

	pitch->fftObj=fftObj;

	pitch->fftLength=fftLength;
	pitch->slideLength=_slideLength;
	pitch->radix2Exp=_radix2Exp;

	pitch->interpFFTLength=interpFFTLength;

	pitch->samplate=_samplate;
	pitch->winType=_winType;

	pitch->lowFre=_lowFre;
	pitch->highFre=_highFre;

	pitch->harmonicCount=_hc;
	pitch->isContinue=_isContinue;

	__pitchLHSObj_initData(pitch);

	return status;
}

int pitchLHSObj_calTimeLength(PitchLHSObj pitchLHSObj,int dataLength){
	int fftLength=0; 
	int slideLength=0;
	int tailDataLength=0;

	int isContinue=0;

	int timeLength=0;

	fftLength=pitchLHSObj->fftLength;
	slideLength=pitchLHSObj->slideLength;
	tailDataLength=pitchLHSObj->tailDataLength;

	isContinue=pitchLHSObj->isContinue;

	if(isContinue){
		dataLength+=tailDataLength; // outTimeLength
	}

	if(dataLength<fftLength){
		return 0;
	}

	timeLength=(dataLength-fftLength)/slideLength+1;
	return timeLength;
}

static void __pitchLHSObj_initData(PitchLHSObj pitchLHSObj){
	int fftLength=0;
	int interpFFTLength=0;

	int minIndex=0; // edge
	int maxIndex=0;

	float *winDataArr=NULL; // fftLength

	int samplate=0;
	WindowType winType=Window_Hamm;

	float lowFre=0;
	float highFre=0;

	int harmonicCount=0;
	int _k=0;

	fftLength=pitchLHSObj->fftLength;
	interpFFTLength=pitchLHSObj->interpFFTLength;

	samplate=pitchLHSObj->samplate;
	winType=pitchLHSObj->winType;

	lowFre=pitchLHSObj->lowFre;
	highFre=pitchLHSObj->highFre;

	harmonicCount=pitchLHSObj->harmonicCount;

	minIndex=ceilf(lowFre);
	maxIndex=floorf(highFre);

	_k=samplate/(maxIndex+1);
	if(harmonicCount>_k){
		harmonicCount=_k;
		if(!harmonicCount){
			harmonicCount=1;
		}
	}

	winDataArr=window_calFFTWindow(winType, fftLength);

	pitchLHSObj->minIndex=minIndex;
	pitchLHSObj->maxIndex=maxIndex;

	pitchLHSObj->harmonicCount=harmonicCount;

	pitchLHSObj->winDataArr=winDataArr;

	pitchLHSObj->realArr1=__vnew(interpFFTLength, NULL);
	pitchLHSObj->imageArr1=__vnew(interpFFTLength, NULL);

	pitchLHSObj->dataArr1=__vnew(interpFFTLength, NULL);
	pitchLHSObj->tailDataArr=__vnew(fftLength, NULL);
}

static int __pitchLHSObj_dealData(PitchLHSObj pitchLHSObj,float *dataArr,int dataLength){
	int status=1;

	int fftLength=0; 
	int slideLength=0;

	int isContinue=0;

	float *tailDataArr=NULL;
	int tailDataLength=0;

	float *curDataArr=NULL;
	int curDataLength=0; 

	int timeLength=0;
	int interpFFTLength=0; 

	int timeLen=0;
	int tailLen=0;

	int totalLength=0;

	fftLength=pitchLHSObj->fftLength;
	slideLength=pitchLHSObj->slideLength;

	isContinue=pitchLHSObj->isContinue;

	tailDataArr=pitchLHSObj->tailDataArr;
	tailDataLength=pitchLHSObj->tailDataLength;

	curDataArr=pitchLHSObj->curDataArr;
	curDataLength=pitchLHSObj->curDataLength;

	timeLength=pitchLHSObj->timeLength;
	interpFFTLength=pitchLHSObj->interpFFTLength;

	if(isContinue){
		totalLength=tailDataLength+dataLength;
	}
	else{
		totalLength=dataLength;
	}

	if(totalLength<fftLength){
		tailLen=totalLength;
		status=0;
	}

	if(status){
		__calTimeAndTailLen(totalLength, fftLength, slideLength, &timeLen, &tailLen);
	}

	if(status){ // has timeLen, cal curDataArr
		if(totalLength>curDataLength||
			curDataLength>2*totalLength){

			free(curDataArr);
			curDataArr=(float *)calloc(totalLength+fftLength, sizeof(float ));
		}

		curDataLength=0;
		if(isContinue&&tailDataLength<0){
			memcpy(curDataArr, dataArr-tailDataLength, (dataLength+tailDataLength)*sizeof(float ));
			curDataLength=(dataLength+tailDataLength);
		}
		else{
			if(isContinue&&tailDataLength>0){ // has & tail
				memcpy(curDataArr, tailDataArr, tailDataLength*sizeof(float ));
				curDataLength+=tailDataLength;
			}

			memcpy(curDataArr+curDataLength, dataArr, dataLength*sizeof(float ));
			curDataLength+=dataLength;
		}

		// tailDataArr
		tailDataLength=0; // reset !!!
		if(isContinue){ 
			if(tailLen>0){
				memcpy(tailDataArr,curDataArr+(curDataLength-tailLen),tailLen*sizeof(float ));
			}
			
			tailDataLength=tailLen;
		}

		// update cache
		if(pitchLHSObj->timeLength<timeLen||
			pitchLHSObj->timeLength>timeLen*2){ 
			free(pitchLHSObj->mDbArr);
			free(pitchLHSObj->mSumArr);
	
			pitchLHSObj->mDbArr=__vnew(timeLen*interpFFTLength,NULL);
			pitchLHSObj->mSumArr=__vnew(timeLen*interpFFTLength,NULL);
		}
	}
	else{
		if(isContinue){ 
			if(tailLen>0){
				if(tailDataLength>=0){
					memcpy(tailDataArr+tailDataLength,dataArr,dataLength*sizeof(float ));
				}
				else{
					memcpy(tailDataArr,dataArr-tailDataLength,(dataLength+tailDataLength)*sizeof(float ));
				}
			}
			
			tailDataLength=tailLen;
		}
		else{
			tailDataLength=0;
		}
	}

	pitchLHSObj->tailDataLength=tailDataLength;

	pitchLHSObj->curDataArr=curDataArr;
	pitchLHSObj->curDataLength=curDataLength;

	pitchLHSObj->timeLength=timeLen;

	return status;
}

void pitchLHSObj_pitch(PitchLHSObj pitchLHSObj,float *dataArr,int dataLength,
					float *freArr){
	int status=0;

	if(!dataArr||dataLength<=0){
		return;
	}

	// 1. deal data
	status=__pitchLHSObj_dealData(pitchLHSObj,dataArr,dataLength);
	if(!status){
		return;
	}

	// 1. calDb
	__pitchLHSObj_calDb(pitchLHSObj);

	// 2. calSum
	__pitchLHSObj_calSum(pitchLHSObj);

	// 3. peak pick
	__pitchLHSObj_dealResult(pitchLHSObj,freArr);

}

static void __pitchLHSObj_calDb(PitchLHSObj pitchLHSObj){
	FFTObj fftObj=NULL;

	int fftLength=0;
	int slideLength=0;
	int timeLength=0;

	float *winDataArr=NULL; // fftLength

	int interpFFTLength=0;

	float *mDbArr=NULL; 

	float *realArr1=NULL; 
	float *imageArr1=NULL;

	float *dataArr1=NULL; 
	float *curDataArr=NULL;

	fftObj=pitchLHSObj->fftObj;

	fftLength=pitchLHSObj->fftLength;
	slideLength=pitchLHSObj->slideLength;
	timeLength=pitchLHSObj->timeLength;

	winDataArr=pitchLHSObj->winDataArr;

	interpFFTLength=pitchLHSObj->interpFFTLength;

	mDbArr=pitchLHSObj->mDbArr;

	realArr1=pitchLHSObj->realArr1;
	imageArr1=pitchLHSObj->imageArr1;

	dataArr1=pitchLHSObj->dataArr1;
	curDataArr=pitchLHSObj->curDataArr;

	for(int i=0;i<timeLength;i++){
		// 0. reset
		memset(dataArr1, 0, sizeof(float )*interpFFTLength);

		// 1. fft
		memcpy(dataArr1, curDataArr+i*slideLength, sizeof(float )*fftLength);
		__vmul(dataArr1, winDataArr, fftLength, NULL);

		fftObj_fft(fftObj, dataArr1, NULL, realArr1, imageArr1);

		__vcabs(realArr1,imageArr1,interpFFTLength,mDbArr+i*interpFFTLength); // __vcsqure
		__vlog(mDbArr+i*interpFFTLength, interpFFTLength, NULL);
	}
}

static void __pitchLHSObj_calSum(PitchLHSObj pitchLHSObj){
	int timeLength=0;
	int interpFFTLength=0;

	int harmonicCount=0;

	int minIndex=0;
	int maxIndex=0;

	float *mDbArr=NULL; 
	float *mSumArr=NULL; 

	float *vArr1=NULL;
	float sum=0;

	timeLength=pitchLHSObj->timeLength;
	interpFFTLength=pitchLHSObj->interpFFTLength;

	harmonicCount=pitchLHSObj->harmonicCount;

	minIndex=pitchLHSObj->minIndex;
	maxIndex=pitchLHSObj->maxIndex;

	mDbArr=pitchLHSObj->mDbArr;
	mSumArr=pitchLHSObj->mSumArr;

	for(int i=0;i<timeLength;i++){
		vArr1=mDbArr+i*interpFFTLength;
		for(int j=0;j<maxIndex+1;j++){
			sum=0;
			for(int k=0;k<harmonicCount;k++){
				sum+=vArr1[j*(k+1)];
			}
			
			mSumArr[i*interpFFTLength+j]=sum;
		}
	}
}

static void __pitchLHSObj_dealResult(PitchLHSObj pitchLHSObj,float *freArr){
	int timeLength=0;
	int interpFFTLength=0;

	int minIndex=0;
	int maxIndex=0;

	int samplate=0;
	float *mSumArr=NULL; 

	float value1=0;
	int index1=0;

	timeLength=pitchLHSObj->timeLength;
	interpFFTLength=pitchLHSObj->interpFFTLength;

	minIndex=pitchLHSObj->minIndex;
	maxIndex=pitchLHSObj->maxIndex;

	samplate=pitchLHSObj->samplate;
	mSumArr=pitchLHSObj->mSumArr;
	for(int i=0;i<timeLength;i++){
		util_peakPick(mSumArr+i*interpFFTLength, maxIndex+1,minIndex, maxIndex, 1, 1, &value1, &index1);
		freArr[i]=(index1+1)*(1.0*samplate/interpFFTLength);
	}
}

void pitchLHSObj_enableDebug(PitchLHSObj pitchLHSObj,int isDebug){

	pitchLHSObj->isDebug=1;
}

void pitchLHSObj_free(PitchLHSObj pitchLHSObj){

	if(pitchLHSObj){
		fftObj_free(pitchLHSObj->fftObj);

		free(pitchLHSObj->winDataArr);

		free(pitchLHSObj->mDbArr);
		free(pitchLHSObj->mSumArr);

		free(pitchLHSObj->realArr1);
		free(pitchLHSObj->imageArr1);

		free(pitchLHSObj->dataArr1);

		free(pitchLHSObj->tailDataArr);
		free(pitchLHSObj->curDataArr);

		free(pitchLHSObj);
	}
}

static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int *timeLength,int *tailLength){
	int timeLen=0;
	int tailLen=0;

	timeLen=(dataLength-fftLength)/slideLength+1;
	tailLen=(dataLength-fftLength)%slideLength+(fftLength-slideLength);

	if(timeLength){
		*timeLength=timeLen;
	}

	if(tailLength){
		*tailLength=tailLen;
	}
}









