// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"

#include "_pitch_cep.h"

struct OpaquePitchCEP{
	int isContinue;

	FFTObj fftObj;

	int fftLength; 
	int slideLength;
	int radix2Exp; // fftLength

	int cepFFTLength; 
	int timeLength;

	int minIndex; // edge
	int maxIndex;

	float *winDataArr; // fftLength
	float *mCepArr; // timeLength*cepFFTLength

	// cache data
	float *realArr1; // cepFFTLength
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

static void __pitchCEPObj_initData(PitchCEPObj pitchCEPObj);
static int __pitchCEPObj_dealData(PitchCEPObj pitchCEPObj,float *dataArr,int dataLength);

static void __pitchCEPObj_calCep(PitchCEPObj pitchCEPObj);
static void __pitchCEPObj_dealResult(PitchCEPObj pitchCEPObj,float *freArr);

/***
	samplate 32000
	lowFre 32,
	highFre 2000
	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int pitchCEPObj_new(PitchCEPObj *pitchCEPObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
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

	FFTObj fftObj=NULL;
	PitchCEPObj pitch=NULL;

	pitch=*pitchCEPObj=(PitchCEPObj )calloc(1,sizeof(struct OpaquePitchCEP ));

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

	fftObj_new(&fftObj, _radix2Exp+1);

	pitch->fftObj=fftObj;

	pitch->fftLength=fftLength;
	pitch->slideLength=_slideLength;
	pitch->radix2Exp=_radix2Exp;

	pitch->cepFFTLength=fftLength*2;

	pitch->samplate=_samplate;
	pitch->winType=_winType;

	pitch->lowFre=_lowFre;
	pitch->highFre=_highFre;

	pitch->isContinue=_isContinue;

	__pitchCEPObj_initData(pitch);

	return status;
}

int pitchCEPObj_calTimeLength(PitchCEPObj pitchCEPObj,int dataLength){
	int fftLength=0; 
	int slideLength=0;
	int tailDataLength=0;

	int isContinue=0;

	int timeLength=0;

	fftLength=pitchCEPObj->fftLength;
	slideLength=pitchCEPObj->slideLength;
	tailDataLength=pitchCEPObj->tailDataLength;

	isContinue=pitchCEPObj->isContinue;

	if(isContinue){
		dataLength+=tailDataLength; // outTimeLength
	}

	if(dataLength<fftLength){
		return 0;
	}

	timeLength=(dataLength-fftLength)/slideLength+1;
	return timeLength;
}

static void __pitchCEPObj_initData(PitchCEPObj pitchCEPObj){
	int fftLength=0;
	int cepFFTLength=0;

	int minIndex=-1; // edge
	int maxIndex=0;

	float *winDataArr=NULL; // fftLength

	int samplate=0;
	WindowType winType=Window_Hamm;

	float lowFre=0;
	float highFre=0;

	fftLength=pitchCEPObj->fftLength;
	cepFFTLength=pitchCEPObj->cepFFTLength;

	samplate=pitchCEPObj->samplate;
	winType=pitchCEPObj->winType;

	lowFre=pitchCEPObj->lowFre;
	highFre=pitchCEPObj->highFre;

	minIndex=roundf(samplate/highFre);
	maxIndex=roundf(samplate/lowFre);

	winDataArr=window_calFFTWindow(winType, fftLength);

	pitchCEPObj->minIndex=minIndex;
	pitchCEPObj->maxIndex=maxIndex;

	pitchCEPObj->winDataArr=winDataArr;

	pitchCEPObj->realArr1=__vnew(cepFFTLength, NULL);
	pitchCEPObj->imageArr1=__vnew(cepFFTLength, NULL);

	pitchCEPObj->realArr2=__vnew(cepFFTLength, NULL);
	pitchCEPObj->imageArr2=__vnew(cepFFTLength, NULL);

	pitchCEPObj->dataArr1=__vnew(cepFFTLength, NULL);
	pitchCEPObj->tailDataArr=__vnew(fftLength, NULL);
}

static int __pitchCEPObj_dealData(PitchCEPObj pitchCEPObj,float *dataArr,int dataLength){
	int status=1;

	int fftLength=0; 
	int slideLength=0;

	int isContinue=0;

	float *tailDataArr=NULL;
	int tailDataLength=0;

	float *curDataArr=NULL;
	int curDataLength=0; 

	int timeLength=0;
	int cepFFTLength=0; 

	int timeLen=0;
	int tailLen=0;

	int totalLength=0;

	fftLength=pitchCEPObj->fftLength;
	slideLength=pitchCEPObj->slideLength;

	isContinue=pitchCEPObj->isContinue;

	tailDataArr=pitchCEPObj->tailDataArr;
	tailDataLength=pitchCEPObj->tailDataLength;

	curDataArr=pitchCEPObj->curDataArr;
	curDataLength=pitchCEPObj->curDataLength;

	timeLength=pitchCEPObj->timeLength;
	cepFFTLength=pitchCEPObj->cepFFTLength;

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
		if(pitchCEPObj->timeLength<timeLen||
			pitchCEPObj->timeLength>timeLen*2){ 
			free(pitchCEPObj->mCepArr);
	
			pitchCEPObj->mCepArr=__vnew(timeLen*cepFFTLength,NULL);
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

	pitchCEPObj->tailDataLength=tailDataLength;

	pitchCEPObj->curDataArr=curDataArr;
	pitchCEPObj->curDataLength=curDataLength;

	pitchCEPObj->timeLength=timeLen;

	return status;
}

void pitchCEPObj_pitch(PitchCEPObj pitchCEPObj,float *dataArr,int dataLength,
					float *freArr){
	int status=0;

	if(!dataArr||dataLength<=0){
		return;
	}

	// 1. deal data
	status=__pitchCEPObj_dealData(pitchCEPObj,dataArr,dataLength);
	if(!status){
		return;
	}

	// 2. cep
	__pitchCEPObj_calCep(pitchCEPObj);

	// 3. peak pick
	__pitchCEPObj_dealResult(pitchCEPObj,freArr);
}

static void __pitchCEPObj_calCep(PitchCEPObj pitchCEPObj){
	FFTObj fftObj=NULL;

	int fftLength=0;
	int slideLength=0;
	int timeLength=0;

	WindowType winType=Window_Rect;
	float *winDataArr=NULL; // fftLength

	int minIndex=0; // edge
	int maxIndex=0;

	int cepFFTLength=0;

	float *mCepArr=NULL; 

	float *realArr1=NULL; 
	float *imageArr1=NULL;

	float *realArr2=NULL; 

	float *dataArr1=NULL; 
	float *curDataArr=NULL;

	fftObj=pitchCEPObj->fftObj;

	fftLength=pitchCEPObj->fftLength;
	slideLength=pitchCEPObj->slideLength;
	timeLength=pitchCEPObj->timeLength;

	winType=pitchCEPObj->winType;
	winDataArr=pitchCEPObj->winDataArr;

	minIndex=pitchCEPObj->minIndex;
	maxIndex=pitchCEPObj->maxIndex;

	cepFFTLength=pitchCEPObj->cepFFTLength;

	mCepArr=pitchCEPObj->mCepArr;

	realArr1=pitchCEPObj->realArr1;
	imageArr1=pitchCEPObj->imageArr1;

	realArr2=pitchCEPObj->realArr2;

	dataArr1=pitchCEPObj->dataArr1;
	curDataArr=pitchCEPObj->curDataArr;

	for(int i=0;i<timeLength;i++){
		// 0. reset
		memset(dataArr1, 0, sizeof(float )*cepFFTLength);

		// 1. cep
		memcpy(dataArr1, curDataArr+i*slideLength, sizeof(float )*fftLength);
		if(winType!=Window_Rect){
			__vmul(dataArr1, winDataArr, fftLength, NULL);
		}
		
		fftObj_fft(fftObj, dataArr1, NULL, realArr1, imageArr1);

		__vcsquare(realArr1,imageArr1,cepFFTLength,realArr2); 
		__vlog(realArr2, cepFFTLength, NULL);

		fftObj_ifft(fftObj, realArr2, NULL, mCepArr+i*cepFFTLength, imageArr1);
	}
}

static void __pitchCEPObj_dealResult(PitchCEPObj pitchCEPObj,float *freArr){
	int timeLength=0;
	int cepFFTLength=0;

	int minIndex=0;
	int maxIndex=0;

	int samplate=0;
	float *mCepArr=NULL; 

	float value1=0;
	int index1=0;

	timeLength=pitchCEPObj->timeLength;
	cepFFTLength=pitchCEPObj->cepFFTLength;

	minIndex=pitchCEPObj->minIndex;
	maxIndex=pitchCEPObj->maxIndex;

	samplate=pitchCEPObj->samplate;
	mCepArr=pitchCEPObj->mCepArr;
	for(int i=0;i<timeLength;i++){
		util_peakPick(mCepArr+i*cepFFTLength, cepFFTLength,minIndex, maxIndex, 1, 1, &value1, &index1);
		freArr[i]=1.0*samplate/(index1+1);
	}
}	

void pitchCEPObj_enableDebug(PitchCEPObj pitchCEPObj,int isDebug){

	pitchCEPObj->isDebug=isDebug;
}

void pitchCEPObj_free(PitchCEPObj pitchCEPObj){

	if(pitchCEPObj){
		fftObj_free(pitchCEPObj->fftObj);

		free(pitchCEPObj->winDataArr);
		free(pitchCEPObj->mCepArr);

		free(pitchCEPObj->realArr1);
		free(pitchCEPObj->imageArr1);

		free(pitchCEPObj->realArr2);
		free(pitchCEPObj->imageArr2);

		free(pitchCEPObj->dataArr1);
		free(pitchCEPObj->tailDataArr);

		free(pitchCEPObj->curDataArr);

		free(pitchCEPObj);
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










