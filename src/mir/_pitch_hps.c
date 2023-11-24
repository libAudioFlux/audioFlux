// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"

#include "_pitch_hps.h"

struct OpaquePitchHPS{
	int isContinue;

	FFTObj fftObj;

	int fftLength; 
	int slideLength;
	int radix2Exp; // fftLength

	int interpFFTLength;
	int timeLength;

	int minIndex; // edge
	int maxIndex;

	int harmonicCount;
	float *winDataArr; // fftLength

	float *mHpsArr; // timeLength*interpFFTLength

	// cache data
	float *realArr1; // interpFFTLength
	float *imageArr1;

	float *realArr2; 
	float *imageArr2;

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

static void __pitchHPSObj_initData(PitchHPSObj pitchHPSObj);
static int __pitchHPSObj_dealData(PitchHPSObj pitchHPSObj,float *dataArr,int dataLength);

static void __pitchHPSObj_calHps(PitchHPSObj pitchHPSObj);
static void __pitchHPSObj_dealResult(PitchHPSObj pitchHPSObj,float *freArr);

/***
	samplate 32000
	lowFre 32,
	highFre 2000
	radix2Exp 12
	WindowType Hamm
	slideLength (1<<radix2Exp)/4
	
	harmonicCount 5 >0
	isContinue 0
****/
int pitchHPSObj_new(PitchHPSObj *pitchHPSObj,
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
	PitchHPSObj pitch=NULL;

	pitch=*pitchHPSObj=(PitchHPSObj )calloc(1,sizeof(struct OpaquePitchHPS ));

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
		if(*windowType<=Window_Hamm){
			_winType=*windowType;
		}
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

	__pitchHPSObj_initData(pitch);

	return status;
}

int pitchHPSObj_calTimeLength(PitchHPSObj pitchHPSObj,int dataLength){
	int fftLength=0; 
	int slideLength=0;
	int tailDataLength=0;

	int isContinue=0;

	int timeLength=0;

	fftLength=pitchHPSObj->fftLength;
	slideLength=pitchHPSObj->slideLength;
	tailDataLength=pitchHPSObj->tailDataLength;

	isContinue=pitchHPSObj->isContinue;

	if(isContinue){
		dataLength+=tailDataLength; // outTimeLength
	}

	if(dataLength<fftLength){
		return 0;
	}

	timeLength=(dataLength-fftLength)/slideLength+1;
	return timeLength;
}

static void __pitchHPSObj_initData(PitchHPSObj pitchHPSObj){
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

	fftLength=pitchHPSObj->fftLength;
	interpFFTLength=pitchHPSObj->interpFFTLength;

	samplate=pitchHPSObj->samplate;
	winType=pitchHPSObj->winType;

	lowFre=pitchHPSObj->lowFre;
	highFre=pitchHPSObj->highFre;

	harmonicCount=pitchHPSObj->harmonicCount;

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

	pitchHPSObj->minIndex=minIndex;
	pitchHPSObj->maxIndex=maxIndex;

	pitchHPSObj->winDataArr=winDataArr;

	pitchHPSObj->realArr1=__vnew(interpFFTLength, NULL);
	pitchHPSObj->imageArr1=__vnew(interpFFTLength, NULL);

	pitchHPSObj->realArr2=__vnew(interpFFTLength, NULL);
	pitchHPSObj->imageArr2=__vnew(interpFFTLength, NULL);

	pitchHPSObj->dataArr1=__vnew(interpFFTLength, NULL);
	pitchHPSObj->tailDataArr=__vnew(fftLength, NULL);
}

static int __pitchHPSObj_dealData(PitchHPSObj pitchHPSObj,float *dataArr,int dataLength){
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

	fftLength=pitchHPSObj->fftLength;
	slideLength=pitchHPSObj->slideLength;

	isContinue=pitchHPSObj->isContinue;

	tailDataArr=pitchHPSObj->tailDataArr;
	tailDataLength=pitchHPSObj->tailDataLength;

	curDataArr=pitchHPSObj->curDataArr;
	curDataLength=pitchHPSObj->curDataLength;

	timeLength=pitchHPSObj->timeLength;
	interpFFTLength=pitchHPSObj->interpFFTLength;

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
		if(pitchHPSObj->timeLength<timeLen||
			pitchHPSObj->timeLength>timeLen*2){ 
			free(pitchHPSObj->mHpsArr);
	
			pitchHPSObj->mHpsArr=__vnew(timeLen*interpFFTLength,NULL);
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

	pitchHPSObj->tailDataLength=tailDataLength;

	pitchHPSObj->curDataArr=curDataArr;
	pitchHPSObj->curDataLength=curDataLength;

	pitchHPSObj->timeLength=timeLen;

	return status;
}


void pitchHPSObj_pitch(PitchHPSObj pitchHPSObj,float *dataArr,int dataLength,
					float *freArr){
	int status=0;

	if(!dataArr||dataLength<=0){
		return;
	}

	// 1. deal data
	status=__pitchHPSObj_dealData(pitchHPSObj,dataArr,dataLength);
	if(!status){
		return;
	}

	// 1. calHps
	__pitchHPSObj_calHps(pitchHPSObj);

	// 2. peak pick
	__pitchHPSObj_dealResult(pitchHPSObj,freArr);

}

static void __pitchHPSObj_calHps(PitchHPSObj pitchHPSObj){
	FFTObj fftObj=NULL;

	int fftLength=0;
	int slideLength=0;
	int timeLength=0;

	WindowType winType=Window_Hamm;
	float *winDataArr=NULL; // fftLength

	int interpFFTLength=0;
	int harmonicCount=0;

	int minIndex=0;
	int maxIndex=0;

	float *mHpsArr=NULL; 
	float hps=0;

	float *realArr1=NULL; 
	float *imageArr1=NULL;

	float *realArr2=NULL; 

	float *dataArr1=NULL; 
	float *curDataArr=NULL;

	fftObj=pitchHPSObj->fftObj;

	fftLength=pitchHPSObj->fftLength;
	slideLength=pitchHPSObj->slideLength;
	timeLength=pitchHPSObj->timeLength;

	winType=pitchHPSObj->winType;
	winDataArr=pitchHPSObj->winDataArr;

	interpFFTLength=pitchHPSObj->interpFFTLength;
	harmonicCount=pitchHPSObj->harmonicCount;

	minIndex=pitchHPSObj->minIndex;
	maxIndex=pitchHPSObj->maxIndex;

	mHpsArr=pitchHPSObj->mHpsArr;

	realArr1=pitchHPSObj->realArr1;
	imageArr1=pitchHPSObj->imageArr1;

	realArr2=pitchHPSObj->realArr2;

	dataArr1=pitchHPSObj->dataArr1;
	curDataArr=pitchHPSObj->curDataArr;

	for(int i=0;i<timeLength;i++){
		// 0. reset
		memset(dataArr1, 0, sizeof(float )*interpFFTLength);

		// 1. fft
		memcpy(dataArr1, curDataArr+i*slideLength, sizeof(float )*fftLength);
		if(winType!=Window_Rect){
			__vmul(dataArr1, winDataArr, fftLength, NULL);
		}
		
		fftObj_fft(fftObj, dataArr1, NULL, realArr1, imageArr1);

		__vcabs(realArr1,imageArr1,interpFFTLength,realArr2); 
		
		// 2.hps
		for(int j=0;j<maxIndex+1;j++){
			hps=1;
			for(int k=0;k<harmonicCount;k++){
				hps*=realArr2[j*(k+1)];
			}

			mHpsArr[i*interpFFTLength+j]=hps;
		}
	}
}

static void __pitchHPSObj_dealResult(PitchHPSObj pitchHPSObj,float *freArr){
	int timeLength=0;
	int interpFFTLength=0;

	int minIndex=0;
	int maxIndex=0;

	int samplate=0;
	float *mHpsArr=NULL; 

	float value1=0;
	int index1=0;

	timeLength=pitchHPSObj->timeLength;
	interpFFTLength=pitchHPSObj->interpFFTLength;

	minIndex=pitchHPSObj->minIndex;
	maxIndex=pitchHPSObj->maxIndex;

	samplate=pitchHPSObj->samplate;
	mHpsArr=pitchHPSObj->mHpsArr;
	for(int i=0;i<timeLength;i++){
		util_peakPick(mHpsArr+i*interpFFTLength, maxIndex+1,minIndex, maxIndex, 1, 1, &value1, &index1);
		freArr[i]=(index1+1)*(1.0*samplate/interpFFTLength);
	}
}

void pitchHPSObj_enableDebug(PitchHPSObj pitchHPSObj,int isDebug){

	pitchHPSObj->isDebug=isDebug;
}

void pitchHPSObj_free(PitchHPSObj pitchHPSObj){

	if(pitchHPSObj){
		fftObj_free(pitchHPSObj->fftObj);

		free(pitchHPSObj->winDataArr);
		free(pitchHPSObj->mHpsArr);

		free(pitchHPSObj->realArr1);
		free(pitchHPSObj->imageArr1);

		free(pitchHPSObj->realArr2);
		free(pitchHPSObj->imageArr2);

		free(pitchHPSObj->dataArr1);
		free(pitchHPSObj->tailDataArr);

		free(pitchHPSObj->curDataArr);

		free(pitchHPSObj);
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













