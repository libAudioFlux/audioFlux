// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"

#include "_pitch_ncf.h"

struct OpaquePitchNCF{
	int isContinue;

	FFTObj fftObj;

	int fftLength; 
	int slideLength;
	int radix2Exp; // fftLength

	int corrFFTLength; 
	int timeLength;

	int minIndex; // edge
	int maxIndex;

	float *winDataArr; // fftLength
	float *mCorrArr; // timeLength*corrFFTLength

	// cache data
	float *realArr1; // corrFFTLength
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

static void __pitchNCFObj_initData(PitchNCFObj pitchNCFObj);
static int __pitchNCFObj_dealData(PitchNCFObj pitchNCFObj,float *dataArr,int dataLength);

static void __pitchNCFObj_calCorr(PitchNCFObj pitchNCFObj);
static void __pitchNCFObj_dealResult(PitchNCFObj pitchNCFObj,float *freArr);

/***
	samplate 32000
	lowFre 32,
	highFre 2000
	radix2Exp 12
	WindowType rect
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int pitchNCFObj_new(PitchNCFObj *pitchNCFObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				int *isContinue){
	int status=0;
	
	float _lowFre=32; // 27.5
	float _highFre=2000; // 2093/4186

	int _samplate=32000;
	int _radix2Exp=12;
	int _slideLength=0;
	WindowType _winType=Window_Rect;
	int _isContinue=0;

	int fftLength=0;

	FFTObj fftObj=NULL;
	PitchNCFObj pitch=NULL;

	pitch=*pitchNCFObj=(PitchNCFObj )calloc(1,sizeof(struct OpaquePitchNCF ));

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

	fftObj_new(&fftObj, _radix2Exp+1);

	pitch->fftObj=fftObj;

	pitch->fftLength=fftLength;
	pitch->slideLength=_slideLength;
	pitch->radix2Exp=_radix2Exp;

	pitch->corrFFTLength=fftLength*2;

	pitch->samplate=_samplate;
	pitch->winType=_winType;

	pitch->lowFre=_lowFre;
	pitch->highFre=_highFre;

	pitch->isContinue=_isContinue;

	__pitchNCFObj_initData(pitch);

	return status;
}

int pitchNCFObj_calTimeLength(PitchNCFObj pitchNCFObj,int dataLength){
	int fftLength=0; 
	int slideLength=0;
	int tailDataLength=0;

	int isContinue=0;

	int timeLength=0;

	fftLength=pitchNCFObj->fftLength;
	slideLength=pitchNCFObj->slideLength;
	tailDataLength=pitchNCFObj->tailDataLength;

	isContinue=pitchNCFObj->isContinue;

	if(isContinue){
		dataLength+=tailDataLength; // outTimeLength
	}

	if(dataLength<fftLength){
		return 0;
	}

	timeLength=(dataLength-fftLength)/slideLength+1;
	return timeLength;
}

static void __pitchNCFObj_initData(PitchNCFObj pitchNCFObj){
	int fftLength=0;
	int corrFFTLength=0;

	int minIndex=-1; // edge
	int maxIndex=0;

	float *winDataArr=NULL; // fftLength

	int samplate=0;
	WindowType winType=Window_Hamm;

	float lowFre=0;
	float highFre=0;

	fftLength=pitchNCFObj->fftLength;
	corrFFTLength=pitchNCFObj->corrFFTLength;

	samplate=pitchNCFObj->samplate;
	winType=pitchNCFObj->winType;

	lowFre=pitchNCFObj->lowFre;
	highFre=pitchNCFObj->highFre;

	minIndex=roundf(samplate/highFre);
	maxIndex=roundf(samplate/lowFre);

	winDataArr=window_calFFTWindow(winType, fftLength);

	pitchNCFObj->minIndex=minIndex;
	pitchNCFObj->maxIndex=maxIndex;

	pitchNCFObj->winDataArr=winDataArr;

	pitchNCFObj->realArr1=__vnew(corrFFTLength, NULL);
	pitchNCFObj->imageArr1=__vnew(corrFFTLength, NULL);

	pitchNCFObj->realArr2=__vnew(corrFFTLength, NULL);
	pitchNCFObj->imageArr2=__vnew(corrFFTLength, NULL);

	pitchNCFObj->dataArr1=__vnew(corrFFTLength, NULL);
	pitchNCFObj->tailDataArr=__vnew(fftLength, NULL);
}

static int __pitchNCFObj_dealData(PitchNCFObj pitchNCFObj,float *dataArr,int dataLength){
	int status=1;

	int fftLength=0; 
	int slideLength=0;

	int isContinue=0;

	float *tailDataArr=NULL;
	int tailDataLength=0;

	float *curDataArr=NULL;
	int curDataLength=0; 

	int timeLength=0;
	int corrFFTLength=0; 

	int timeLen=0;
	int tailLen=0;

	int totalLength=0;

	fftLength=pitchNCFObj->fftLength;
	slideLength=pitchNCFObj->slideLength;

	isContinue=pitchNCFObj->isContinue;

	tailDataArr=pitchNCFObj->tailDataArr;
	tailDataLength=pitchNCFObj->tailDataLength;

	curDataArr=pitchNCFObj->curDataArr;
	curDataLength=pitchNCFObj->curDataLength;

	timeLength=pitchNCFObj->timeLength;
	corrFFTLength=pitchNCFObj->corrFFTLength;

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
		if(pitchNCFObj->timeLength<timeLen||
			pitchNCFObj->timeLength>timeLen*2){ 
			free(pitchNCFObj->mCorrArr);
	
			pitchNCFObj->mCorrArr=__vnew(timeLen*corrFFTLength,NULL);
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

	pitchNCFObj->tailDataLength=tailDataLength;

	pitchNCFObj->curDataArr=curDataArr;
	pitchNCFObj->curDataLength=curDataLength;

	pitchNCFObj->timeLength=timeLen;

	return status;
}

void pitchNCFObj_pitch(PitchNCFObj pitchNCFObj,float *dataArr,int dataLength,
					float *freArr){
	int status=0;

	if(!dataArr||dataLength<=0){
		return;
	}

	// 1. deal data
	status=__pitchNCFObj_dealData(pitchNCFObj,dataArr,dataLength);
	if(!status){
		return;
	}

	// 2. corr
	__pitchNCFObj_calCorr(pitchNCFObj);

	// 3. peak pick
	__pitchNCFObj_dealResult(pitchNCFObj,freArr);

}

static void __pitchNCFObj_calCorr(PitchNCFObj pitchNCFObj){
	FFTObj fftObj=NULL;

	int fftLength=0;
	int slideLength=0;
	int timeLength=0;

	WindowType winType=Window_Rect;
	float *winDataArr=NULL; // fftLength

	int minIndex=0; // edge
	int maxIndex=0;

	int corrFFTLength=0;

	float *mCorrArr=NULL; 

	float *realArr1=NULL; 
	float *imageArr1=NULL;

	float *realArr2=NULL; 

	float *dataArr1=NULL; 
	float *curDataArr=NULL;

	int len=0;

	int lagNum=0;
	int padNum=0;

	fftObj=pitchNCFObj->fftObj;

	fftLength=pitchNCFObj->fftLength;
	slideLength=pitchNCFObj->slideLength;
	timeLength=pitchNCFObj->timeLength;

	winType=pitchNCFObj->winType;
	winDataArr=pitchNCFObj->winDataArr;

	minIndex=pitchNCFObj->minIndex;
	maxIndex=pitchNCFObj->maxIndex;

	corrFFTLength=pitchNCFObj->corrFFTLength;

	mCorrArr=pitchNCFObj->mCorrArr;

	realArr1=pitchNCFObj->realArr1;
	imageArr1=pitchNCFObj->imageArr1;

	realArr2=pitchNCFObj->realArr2;

	dataArr1=pitchNCFObj->dataArr1;
	curDataArr=pitchNCFObj->curDataArr;

	len=(maxIndex<corrFFTLength-1?maxIndex:corrFFTLength-1);
	lagNum=(2*len+1)-(minIndex+maxIndex);
	padNum=minIndex-1;
	for(int i=0;i<timeLength;i++){
		float rms=0;

		// 0. reset
		memset(dataArr1, 0, sizeof(float )*corrFFTLength);

		// 1. corr
		memcpy(dataArr1, curDataArr+i*slideLength, sizeof(float )*fftLength);
		if(winType!=Window_Rect){
			__vmul(dataArr1, winDataArr, fftLength, NULL);
		}
		
		fftObj_fft(fftObj, dataArr1, NULL, realArr1, imageArr1);

		__vcsquare(realArr1,imageArr1,corrFFTLength,realArr2); 
		fftObj_ifft(fftObj, realArr2, NULL, realArr1, imageArr1);

		__vmul_value(realArr1, 1.0/sqrtf(corrFFTLength), corrFFTLength, NULL);

		memcpy(realArr2,realArr1+(corrFFTLength-len),sizeof(float )*len);
		memcpy(realArr2+len,realArr1,sizeof(float )*(len+1));

		// 2. norm
		rms=sqrtf(realArr2[maxIndex]);

		memset(mCorrArr+i*corrFFTLength, 0, sizeof(float )*padNum);
		memcpy(mCorrArr+(i*corrFFTLength+padNum),realArr2+(minIndex+maxIndex),sizeof(float )*lagNum);

		__vmul_value(mCorrArr+(i*corrFFTLength+padNum), 1.0/rms, lagNum, NULL);
	}
}

static void __pitchNCFObj_dealResult(PitchNCFObj pitchNCFObj,float *freArr){
	int timeLength=0;
	int corrFFTLength=0;

	int minIndex=0;
	int maxIndex=0;

	int samplate=0;
	float *mCorrArr=NULL; 

	float value1=0;
	int index1=0;

	timeLength=pitchNCFObj->timeLength;
	corrFFTLength=pitchNCFObj->corrFFTLength;

	minIndex=pitchNCFObj->minIndex;
	maxIndex=pitchNCFObj->maxIndex;

	samplate=pitchNCFObj->samplate;
	mCorrArr=pitchNCFObj->mCorrArr;
	for(int i=0;i<timeLength;i++){
		util_peakPick(mCorrArr+i*corrFFTLength, maxIndex+1,minIndex, maxIndex, 1, 1, &value1, &index1);
		freArr[i]=1.0*samplate/(index1+1);
	}
}

void pitchNCFObj_enableDebug(PitchNCFObj pitchNCFObj,int isDebug){

	pitchNCFObj->isDebug=isDebug;
}

void pitchNCFObj_free(PitchNCFObj pitchNCFObj){

	if(pitchNCFObj){
		fftObj_free(pitchNCFObj->fftObj);

		free(pitchNCFObj->winDataArr);
		free(pitchNCFObj->mCorrArr);

		free(pitchNCFObj->realArr1);
		free(pitchNCFObj->imageArr1);

		free(pitchNCFObj->realArr2);
		free(pitchNCFObj->imageArr2);

		free(pitchNCFObj->dataArr1);
		free(pitchNCFObj->tailDataArr);

		free(pitchNCFObj->curDataArr);

		free(pitchNCFObj);
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









