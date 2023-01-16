// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/flux_correct.h"

#include "filterbank/auditory_filterBank.h"

#include "temporal_algorithm.h"
#include "reassign_algorithm.h"

#include "bft_algorithm.h"

struct OpaqueBFT{
	ReassignObj reassignObj;

	int fftLength; // fftLength,timeLength,num
	int timeLength;

	int num;
	float *mFilterBankArr; // num*(fftLength/2+1)
	float *mArr1; // image cache 

	float *freBandArr;
	int *binBandArr;

	TemporalObj tempObj;

	// rea
	float *mRealArr; // rea r,i(timeLength*(fftLength/2+1)) -> power r(timeLength*(fftLength/2+1)) 
	float *mImageArr;

	// params
	int samplate; // filterBank 
	float lowFre;
	float highFre;

	int lowIndex; // only for linear
	int highIndex;

	int binPerOctave; // log

	int radix2Exp; // rea
	WindowType windowType;
	int slideLength;

	SpectralDataType dataType; // data&&filterBank
	SpectralFilterBankScaleType filterScaleType;
	SpectralFilterBankStyleType filterStyleType;
	SpectralFilterBankNormalType filterNormalType;

	int resultType; // 0 complex 1 real->mRealArr
	float normValue;  // default 1; >0 mag=>result power=>S

	int isReassign; 
	int isTemporal; 

};

static void __bftObj_init(BFTObj bftObj);

/***
	num>=2&&num<=2048
	samplate 32000
	lowFre linear 0 mel/bark/erb/log 27.5
	highFre 16000
	binPerOctave 12 >=4&&<=48
	
	radix2Exp 12 
	WindowType "hann"
	slideLength 1024

	filterScaleType "linear"
	filterStyleType "slaney"
	filterNormalType "none"
	dataType "power"

	isReassign 0
****/
int bftObj_new(BFTObj *bftObj,int num,int radix2Exp,
			int *samplate,float *lowFre,float *highFre,int *binPerOctave,
			WindowType *windowType,int *slideLength,
			SpectralFilterBankScaleType *filterScaleType,
			SpectralFilterBankStyleType *filterStyleType,
			SpectralFilterBankNormalType *filterNormalType,
			SpectralDataType *dataType,
			int *isReassign,
			int *isTemporal){
	int status=0;
	BFTObj bft=NULL;

	int fftLength=0;

	int _samplate=32000; // filterBank 
	float _lowFre=0;
	float _highFre=0;

	int lowIndex=0; // only for linear
	int highIndex=0;

	int _binPerOctave=12; // log
	int _radix2Exp=12; // reassign

	WindowType _windowType=Window_Hann;
	int _slideLength=0;

	SpectralDataType _dataType=SpectralData_Power; // data&&filterBank
	SpectralFilterBankScaleType _filterScaleType=SpectralFilterBankScale_Linear;
	SpectralFilterBankStyleType _filterStyleType=SpectralFilterBankStyle_Slaney;
	SpectralFilterBankNormalType _filterNormalType=SpectralFilterBankNormal_None;

	int _isReassign=0;
	int _isTemporal=0;

	if(radix2Exp){
		_radix2Exp=radix2Exp;
		if(_radix2Exp<1||_radix2Exp>30){
			status=-100;
			printf("radix2Exp is error!\n");
			return status;
		}
	}

	fftLength=1<<_radix2Exp;
	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
		}
	}

	if(dataType){
		_dataType=*dataType;
	}
	
	if(filterScaleType){
		_filterScaleType=*filterScaleType;
		if(_filterScaleType>SpectralFilterBankScale_Log){
			printf("scaleType is error!\n");
			return 1;
		}
	}

	if(filterStyleType){
		_filterStyleType=*filterStyleType;
	}
	
	if(filterNormalType){
		_filterNormalType=*filterNormalType;
	}

	_highFre=_samplate/2.0;

	if(lowFre){
		if(*lowFre>=0&&*lowFre<_samplate/2.0){
			_lowFre=*lowFre;
		}
	}

	if(_lowFre==0){
		if(_filterScaleType==SpectralFilterBankScale_Octave||
			_filterScaleType==SpectralFilterBankScale_Log){ // Log/Logspace

			_lowFre=powf(2, -45/12.0)*440;
			_highFre=powf(2, 38/12.0)*440;
		}
	}

	if(highFre){
		if(*highFre>0&&*highFre<=_samplate/2.0){
			_highFre=*highFre;
		}
	}

	if(_highFre<_lowFre){
		_lowFre=0;
		_highFre=_samplate/2.0;
		if(_filterScaleType==SpectralFilterBankScale_Octave||
			_filterScaleType==SpectralFilterBankScale_Log){ // Log/Logspace

			_lowFre=powf(2, -45/12.0)*440;
			_highFre=powf(2, 38/12.0)*440;
		}
	}

	if(binPerOctave){
		if(*binPerOctave>=4&&*binPerOctave<=48){
			_binPerOctave=*binPerOctave;
		}
	}

	if(windowType){
		_windowType=*windowType;
	}

	_slideLength=fftLength/4;
	if(slideLength){
		if(*slideLength>0){ // &&*slideLength<=fftLength ->support not overlap
			_slideLength=*slideLength;
		}
	}

	if(_filterScaleType==SpectralFilterBankScale_Linear){ // linear
		float det=0;

		det=_samplate/(float )fftLength;
		auditory_reviseLinearFre(num, _lowFre,_highFre, det,1,&_lowFre, &_highFre);

		lowIndex=roundf(_lowFre/det);
		highIndex=roundf(_highFre/det);

		if(_highFre>_samplate/2.0){
			printf("scale linear: lowFre and num is large, overflow error\n");
			return -1;
		}
	}
	else if(_filterScaleType==SpectralFilterBankScale_Octave){ // log
		auditory_reviseLogFre(num, _lowFre,_highFre, _binPerOctave,1, &_lowFre, &_highFre);

		if(_highFre>_samplate/2.0){
			printf("scale log: lowFre and num is large, overflow error!\n");
			return -1;
		}
	}

	if(isReassign){
		_isReassign=*isReassign;
	}

	if(isTemporal){
		_isTemporal=*isTemporal;
	}

	if(num<2||num>fftLength/2+1){ 
		printf("num is error!\n");
		return -1;
	}

	bft=*bftObj=(BFTObj )calloc(1, sizeof(struct OpaqueBFT ));

	bft->fftLength=fftLength;
	bft->num=num;

	bft->samplate=_samplate;
	bft->lowFre=_lowFre;
	bft->highFre=_highFre;

	bft->lowIndex=lowIndex;
	bft->highIndex=highIndex;

	bft->binPerOctave=_binPerOctave;

	bft->radix2Exp=_radix2Exp;
	bft->windowType=_windowType;
	bft->slideLength=_slideLength;

	bft->dataType=_dataType;
	bft->filterScaleType=_filterScaleType;
	bft->filterStyleType=_filterStyleType;
	bft->filterNormalType=_filterNormalType;

	bft->normValue=1;

	bft->isReassign=_isReassign;
	bft->isTemporal=_isTemporal;

	__bftObj_init(bft);

	return status;
}

static void __bftObj_init(BFTObj bftObj){
	ReassignObj reassignObj=NULL;

	int num=0;
	int fftLength=0;

	float *mFilterBankArr=NULL;
	float *mArr1=NULL;

	float *freBandArr=NULL;
	int *binBandArr=NULL;

	int samplate=0; // filterBank 

	float lowFre=0;
	float highFre=0;

	int binPerOctave=0;
	int radix2Exp=0; // stft
	WindowType windowType;
	int slideLength=0;

	SpectralFilterBankScaleType filterScaleType;
	SpectralFilterBankStyleType filterStyleType;
	SpectralFilterBankNormalType filterNormalType;

	int isReassign=0;

	TemporalObj tempObj=NULL;

	ReassignType reType=Reassign_None;

	num=bftObj->num;
	fftLength=bftObj->fftLength;

	samplate=bftObj->samplate;
	lowFre=bftObj->lowFre;
	highFre=bftObj->highFre;

	int lowIndex=0; // only for linear
	int highIndex=0;

	binPerOctave=bftObj->binPerOctave;

	radix2Exp=bftObj->radix2Exp;
	windowType=bftObj->windowType;
	slideLength=bftObj->slideLength;

	filterScaleType=bftObj->filterScaleType;
	filterStyleType=bftObj->filterStyleType;
	filterNormalType=bftObj->filterNormalType;

	isReassign=bftObj->isReassign;

	if(isReassign){
		reType=Reassign_All;
	}

	// 1. reassign
	reassignObj_new(&reassignObj,radix2Exp,
					&samplate,&windowType,&slideLength,
					&reType,NULL,
					NULL,NULL);

	// 2. filterBank
	if(filterScaleType>=SpectralFilterBankScale_Linear&&
		filterScaleType<=SpectralFilterBankScale_Log){

		mFilterBankArr=__vnew(num*(fftLength/2+1), NULL);
		mArr1=__vnew(num*(fftLength/2+1), NULL);

		freBandArr=__vnew(num+2, NULL); 
		binBandArr=__vnewi(num+2, NULL);

		if(filterScaleType==SpectralFilterBankScale_Linear){
			float det=0;

			det=samplate/(float )fftLength;

			lowIndex=bftObj->lowIndex;
			highIndex=bftObj->highIndex;

			for(int i=lowIndex,j=0;i<=highIndex;i++,j++){
				freBandArr[j]=i*det;
				binBandArr[j]=i;
			}
		}
		else{
			auditory_filterBank(num,fftLength,samplate,0,
							filterScaleType,filterStyleType,filterNormalType,
							lowFre,highFre,binPerOctave,
							mFilterBankArr,
							freBandArr,
							binBandArr);
		}
	}

	// 3. temproal
	if(bftObj->isTemporal){
		temporalObj_new(&tempObj,&fftLength,&slideLength,&windowType);
	}

	bftObj->reassignObj=reassignObj;

	bftObj->mFilterBankArr=mFilterBankArr;
	bftObj->mArr1=mArr1;

	bftObj->freBandArr=freBandArr;
	bftObj->binBandArr=binBandArr;

	bftObj->tempObj=tempObj;
}

/***
	1. S real/complex
	2. S mag/power/p
	S.dot(filterBank) t*(fftLength/2+1) @ num*(fftLength/2+1) => t*num
	
****/
void bftObj_bft(BFTObj bftObj,float *dataArr,int dataLength,float *mRealArr3,float *mImageArr3){
	ReassignObj reassignObj=NULL;

	int fftLength=0; // fftLength,timeLength,num
	int timeLength=0;

	int num=0;
	float *mFilterBankArr=NULL;

	SpectralDataType dataType;
	SpectralFilterBankScaleType filterScaleType;

	float normValue=1;
	int resultType=0; // 0 complex 1 real->mRealArr

	float *mRealArr=NULL;
	float *mImageArr=NULL;

	int lowIndex=0; // only for linear
	int highIndex=0;

	float *mArr1=NULL;

	reassignObj=bftObj->reassignObj;

	filterScaleType=bftObj->filterScaleType;
	dataType=bftObj->dataType;

	normValue=bftObj->normValue;
	resultType=bftObj->resultType;

	fftLength=bftObj->fftLength;
	
	num=bftObj->num;

	mFilterBankArr=bftObj->mFilterBankArr; // num*(fftLength/2+1)
	mArr1=bftObj->mArr1;

	mRealArr=bftObj->mRealArr; // timeLength*fftLength
	mImageArr=bftObj->mImageArr;

	lowIndex=bftObj->lowIndex; 
	highIndex=bftObj->highIndex;

	timeLength=reassignObj_calTimeLength(reassignObj,dataLength);
	if(bftObj->timeLength<timeLength||
		bftObj->timeLength>timeLength*2){ // 更新缓存
		free(mRealArr);
		free(mImageArr);
		
		mRealArr=__vnew(timeLength*fftLength, NULL);
		mImageArr=__vnew(timeLength*fftLength, NULL);
	}

	// 1. reassign -> timeLength*(fftLength/2+1)
	reassignObj_reassign(reassignObj,dataArr,dataLength,
						mRealArr,mImageArr,
						NULL,NULL);

	// 2. dot
	if(!resultType){ // complex
		if(dataType==SpectralData_Power){
			for(int i=0;i<timeLength*(fftLength/2+1);i++){
				float r1=0;
				float i1=0;

				r1=mRealArr[i];
				i1=mImageArr[i];

				mRealArr[i]=r1*r1-i1*i1;
				mImageArr[i]=2*r1*i1;
			}
		}

		// dot
		if(filterScaleType==SpectralFilterBankScale_Linear){
			for(int i=0;i<timeLength;i++){
				for(int j=lowIndex,k=0;j<=highIndex;j++,k++){
					mRealArr3[i*num+k]=mRealArr[i*(fftLength/2+1)+j];
					mImageArr3[i*num+k]=mImageArr[i*(fftLength/2+1)+j];
				}
			}
		}
		else{
			__mcdot1(mRealArr,mImageArr,
					mFilterBankArr,mArr1,
					timeLength,fftLength/2+1,
					num,fftLength/2+1,
					mRealArr3,mImageArr3);
		}
	}
	else{
		__mcsquare(mRealArr, mImageArr, timeLength, fftLength/2+1,1, mRealArr); // S^2
		
		// S/S^2 -> norm
		if(dataType==SpectralData_Mag){

			for(int i=0;i<timeLength*(fftLength/2+1);i++){
				mRealArr[i]=sqrtf(mRealArr[i]);
			}
		}
		else if(dataType==SpectralData_Power){
			if(normValue!=1){ 
				for(int i=0;i<timeLength*(fftLength/2+1);i++){
					mRealArr[i]=powf(mRealArr[i], normValue);
				}
			}
		}

		// dot
		if(filterScaleType==SpectralFilterBankScale_Linear){
			for(int i=0;i<timeLength;i++){
				for(int j=lowIndex,k=0;j<=highIndex;j++,k++){
					mRealArr3[i*num+k]=mRealArr[i*(fftLength/2+1)+j];
				}
			}
		}
		else{
			__mdot1(mRealArr,mFilterBankArr,
					timeLength,fftLength/2+1,
					num,fftLength/2+1,
					mRealArr3);
		}
		
		// norm
		if(dataType==SpectralData_Mag){
			if(normValue!=1){
				for(int i=0;i<timeLength*num;i++){
					mRealArr3[i]=powf(mRealArr3[i],normValue);
				}
			}
		}
	}

	// 3. temporal
	if(bftObj->isTemporal){
		temporalObj_temporal(bftObj->tempObj,dataArr,dataLength);
	}

	bftObj->timeLength=timeLength;

	bftObj->mRealArr=mRealArr;
	bftObj->mImageArr=mImageArr;
}

// energy/rms/zeroCrossRate
void bftObj_getTemporalData(BFTObj bftObj,float **eArr,float **rArr,float **zArr){

	if(bftObj->tempObj){
		temporalObj_getData(bftObj->tempObj,eArr,rArr,zArr,NULL);
	}
}

int bftObj_calTimeLength(BFTObj bftObj,int dataLength){
	int timeLength=0;

	timeLength=reassignObj_calTimeLength(bftObj->reassignObj,dataLength);
	return timeLength;
}

float *bftObj_getFreBandArr(BFTObj bftObj){

	return bftObj->freBandArr;
}

int *bftObj_getBinBandArr(BFTObj bftObj){

	return bftObj->binBandArr;
}

// 0 complex,1 real ->mRealArr
void bftObj_setResultType(BFTObj bftObj,int type){

	bftObj->resultType=type;
}

void bftObj_setDataNormValue(BFTObj bftObj,float normValue){

	if(normValue>0){
		bftObj->normValue=normValue;
	}
}

void bftObj_free(BFTObj bftObj){
	ReassignObj reassignObj=NULL;

	float *mFilterBankArr=NULL; // num*(fftLength/2+1)

	float *freBandArr=NULL;
	int *binBandArr=NULL;

	float *mRealArr=NULL; // rea r,i(timeLength*(fftLength/2+1)) -> power r(timeLength*(fftLength/2+1)) 
	float *mImageArr=NULL;

	float *mArr1=NULL; 

	TemporalObj tempObj=NULL;

	if(bftObj){
		reassignObj=bftObj->reassignObj;

		mFilterBankArr=bftObj->mFilterBankArr;

		freBandArr=bftObj->freBandArr;
		binBandArr=bftObj->binBandArr;

		mRealArr=bftObj->mRealArr;
		mImageArr=bftObj->mImageArr;

		mArr1=bftObj->mArr1;

		tempObj=bftObj->tempObj;

		reassignObj_free(reassignObj);

		free(mFilterBankArr);

		free(freBandArr);
		free(binBandArr);

		free(mRealArr);
		free(mImageArr);

		free(mArr1);

		temporal_free(tempObj);

		free(bftObj);
	}
}









