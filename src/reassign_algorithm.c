// clang 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/flux_window.h"

#include "stft_algorithm.h"
#include "reassign_algorithm.h"

struct OpaqueReassign{
	int isContinue; 

	ReassignType resType;
	int isPadding;
	
	int samplate;
	STFTObj stftObj; // 基于STFT linear的reassign ???
	
	int slideLength;
	int fftLength; 
	int timeLength;

	float *freArr;  // fftLength/2+1
	float *timeArr; // timeLength

	float *mReFreArr; // timeLength*(fftLength/2+1)
	float *mReTimeArr;

	// win数据相关 fftLength
	float *winArr; // S_h
	float *winDerivativeArr; // S_dh +1
	float *winWeightArr; // s_th

	// stft相关 timeLength*fftLength => timeLength*(fftLength/2+1)
	float *mRealArr1; // S_h
	float *mImageArr1;

	float *mRealArr2; // S_dh
	float *mImageArr2;

	float *mRealArr3; // S_th
	float *mImageArr3;

	// rearrage相关 timeLength*(fftLength/2+1)
	int *mTimeIndexArr;
	int *mFreIndexArr;

	int *mTempIndexArr;

	float thresh; // 0.001

	int resultType; // 0 complex 1 mRealArr1 -> amp
	int order; // >=1 

};

static void _reassignObj_initWindowData(ReassignObj reassignObj);

static void _reassignObj_dealTCorrData(ReassignObj reassignObj,int dataLength);

static void _reassignObj_stft(ReassignObj reassignObj,float *dataArr,int dataLength);
static void _reassignObj_reassignTimeFre(ReassignObj reassignObj,float *dataArr,int dataLength);
static void _reassignObj_filterTimeFre(ReassignObj reassignObj,int dataLength);
static void _reassignObj_rearrage(ReassignObj reassignObj,float *mRealArr4,float *mImageArr4,float *mRealArr5,float *mImageArr5);

/***
	radix2Exp 12
	samplate 32000
	windowType 'hann'
	slideLength fftLength/4

	reType Reassign_All
	thresh 0.001 

	isPadding 0
	isContinue 0
****/
int reassignObj_new(ReassignObj *reassignObj,int radix2Exp,
					int *samplate,WindowType *windowType,int *slideLength,
					ReassignType *reType,float *thresh,
					int *isPadding,int *isContinue){
	int status=0;
	ReassignObj re=NULL;

	ReassignType _resType=Reassign_All;
	int _samplate=32000;
	int _isPadding=0;

	int _radix2Exp=12;
	WindowType _windowType=Window_Hann;
	int _slideLength=0;
	int _isContinue=0;

	STFTObj stftObj=NULL;

	int fftLength=0;

	float *freArr=NULL;  // fftLength/2+1

	float _thresh=0.001;

	re=*reassignObj=(ReassignObj )calloc(1, sizeof(struct OpaqueReassign ));

	if(reType){
		_resType=*reType;
	}

	if(samplate){
		if(*samplate>0){
			_samplate=*samplate;
		}
	}

	if(isPadding){
		_isPadding=*isPadding;
	}

	if(radix2Exp>1&&radix2Exp<31){
		_radix2Exp=radix2Exp;
	}

	if(windowType){
		_windowType=*windowType;
	}

	fftLength=(1<<_radix2Exp);
	if(slideLength){
		if(*slideLength>0){
			_slideLength=*slideLength;
		}
	}

	if(!_slideLength){
		_slideLength=fftLength/4;
	}

	if(thresh){
		if(*thresh>=0){
			_thresh=*thresh;
		}
	}

	stftObj_new(&stftObj, _radix2Exp, &_windowType, &_slideLength, &_isContinue);
	stftObj_enablePadding(stftObj, _isPadding);

	re->isContinue=_isContinue;

	re->resType=_resType;
	re->isPadding=_isPadding;

	re->samplate=_samplate;
	re->stftObj=stftObj;

	re->slideLength=_slideLength;
	re->fftLength=fftLength;

	re->thresh=_thresh;

	// init winArr相关
	_reassignObj_initWindowData(re);

	// fre相关
	freArr=__vlinspace(0, _samplate/2.0, fftLength/2+1, 0);
	re->freArr=freArr;

	return status;
}

int reassignObj_calTimeLength(ReassignObj reassignObj,int dataLength){

	return stftObj_calTimeLength(reassignObj->stftObj, dataLength);
}

// 0 complex,1 mRealArr1 -> amp
void reassignObj_setResultType(ReassignObj reassignObj,int type){

	reassignObj->resultType=type;
}

// order >=1
void reassignObj_setOrder(ReassignObj reassignObj,int order){

	reassignObj->order=order;
}

/***
	1. tcorr
	2. stft => S_h
	3. stft => S_dh&S_th
	4. filter tcorr/fcorr 
	5. rearrage 
****/
void reassignObj_reassign(ReassignObj reassignObj,float *dataArr,int dataLength,
						float *mRealArr4,float *mImageArr4,
						float *mRealArr5,float *mImageArr5){
	if(reassignObj->resType!=Reassign_None){
		// 1. tcorr
		_reassignObj_dealTCorrData(reassignObj, dataLength);

		// 2. S_h
		_reassignObj_stft(reassignObj, dataArr, dataLength);

		// 3. S_dh&S_th
		_reassignObj_reassignTimeFre(reassignObj,dataArr,dataLength);

		// 4. filter
		_reassignObj_filterTimeFre(reassignObj,dataLength);

		// 5. rearrage
		_reassignObj_rearrage(reassignObj,mRealArr4,mImageArr4,mRealArr5,mImageArr5);
	}
	else{
		int fftLength=0; 
		int timeLength=0;

		float *mRealArr1=NULL; // S_h
		float *mImageArr1=NULL;

		fftLength=reassignObj->fftLength;
		timeLength=reassignObj->timeLength;

		mRealArr1=reassignObj->mRealArr1;
		mImageArr1=reassignObj->mImageArr1;

		timeLength=stftObj_calTimeLength(reassignObj->stftObj,dataLength);
		if(reassignObj->timeLength<timeLength||
			reassignObj->timeLength>timeLength*2){
			free(mRealArr1);
			free(mImageArr1);

			mRealArr1=__vnew(timeLength*fftLength, NULL);
			mImageArr1=__vnew(timeLength*fftLength, NULL);

			reassignObj->timeLength=timeLength;
			reassignObj->mRealArr1=mRealArr1;
			reassignObj->mImageArr1=mImageArr1;
		}

		_reassignObj_stft(reassignObj, dataArr, dataLength);

		memcpy(mRealArr4, mRealArr1, sizeof(float )*timeLength*(fftLength/2+1));
		memcpy(mImageArr4, mImageArr1, sizeof(float )*timeLength*(fftLength/2+1));
	}

}

/***
	fre*timeLength -->timeLength*fre
	rearrage use mag^2/(∑win^2) or mag/power???
****/
static void _reassignObj_rearrage(ReassignObj reassignObj,float *mRealArr4,float *mImageArr4,float *mRealArr5,float *mImageArr5){
	int fftLength=0; 
	int timeLength=0;

	float *freArr=NULL;  // fftLength/2+1
	float *timeArr=NULL;

	float *mReFreArr=NULL; // timeLength*(fftLength/2+1)
	float *mReTimeArr=NULL;

	float *mRealArr1=NULL; // S_h
	float *mImageArr1=NULL;

	int *mTimeIndexArr=NULL; 
	int *mFreIndexArr=NULL;

	int *mTempIndexArr=NULL;

	int order=1;

	int totalLength=0; // timeLength*(fftLength/2+1);

	float fmin=0;
	float fmax=0;

	float tmin=0;
	float tmax=0;

	fftLength=reassignObj->fftLength;
	timeLength=reassignObj->timeLength;

	freArr=reassignObj->freArr;
	timeArr=reassignObj->timeArr;

	mReFreArr=reassignObj->mReFreArr;
	mReTimeArr=reassignObj->mReTimeArr;

	mRealArr1=reassignObj->mRealArr1;
	mImageArr1=reassignObj->mImageArr1;

	mTimeIndexArr=reassignObj->mTimeIndexArr;
	mFreIndexArr=reassignObj->mFreIndexArr;

	mTempIndexArr=reassignObj->mTempIndexArr;

	totalLength=timeLength*(fftLength/2+1);

	order=reassignObj->order;

	fmin=freArr[0];
	fmax=freArr[fftLength/2];

	tmin=timeArr[0];
	tmax=timeArr[timeLength-1];

	// memset(rowIndexArr, 0, sizeof(int )*totalLength);
	if(timeLength>1){
		for(int i=0;i<totalLength;i++){
			mTimeIndexArr[i]=roundf((mReTimeArr[i]-tmin)*(timeLength-1)/(tmax-tmin));
		}
	}

	for(int i=0;i<totalLength;i++){
		mFreIndexArr[i]=roundf((mReFreArr[i]-fmin)*(fftLength/2)/(fmax-fmin));
	}

	// debug
	{
		// printf("mReFreArr is:\n");
		// __mdebug(mReFreArr, timeLength, fftLength/2+1, 1);
		// printf("\n");

		// printf("mReTimeArr is:\n");
		// __mdebug(mReTimeArr, timeLength, fftLength/2+1, 1);
		// printf("\n");

		// printf("timeArr is :\n");
		// __vdebug(timeArr, timeLength, 1);
		// printf("\n");

		// printf("freArr is :\n");
		// __vdebug(freArr, fftLength/2+1, 1);
		// printf("\n");

		// printf("mFreIndexArr is:\n");
		// __mdebugi(mFreIndexArr, timeLength, fftLength/2+1, 1);
		// printf("\n");
	}

	// order
	if(order>1){
		int v=0;

		memset(mTempIndexArr,0,sizeof(int )*totalLength);
		for(int k=0;k<reassignObj->order-1;k++){

			for(int i=0;i<timeLength;i++){
				for(int j=0;j<fftLength/2+1;j++){
					
					v=mFreIndexArr[i*(fftLength/2+1)+j];
					if(v>=0&&v<fftLength/2+1){
						mTempIndexArr[i*(fftLength/2+1)+j]=mFreIndexArr[i*(fftLength/2+1)+v];
					}
				}
			}

			memcpy(mFreIndexArr, mTempIndexArr, sizeof(int )*timeLength*(fftLength/2+1));
			// memset(mTempIndexArr,0,sizeof(int )*timeLength*(fftLength/2+1));
		}
	}
	
	/***
		1. abs/power,df&dt
		2. complex,df&timeArr,modified stft,  
	****/
	for(int i=0;i<timeLength;i++){
		int i1=0;
		int j1=0;

		float v1=0;
		float v2=0;

		for(int j=0;j<fftLength/2+1;j++){
			i1=mTimeIndexArr[i*(fftLength/2+1)+j];
			j1=mFreIndexArr[i*(fftLength/2+1)+j];

			v1=mRealArr1[i*(fftLength/2+1)+j];
			v2=mImageArr1[i*(fftLength/2+1)+j];

			if(j%2==1){
				v1=-v1;
				v2=-v2;
			}

			if(i1>=0&&i1<timeLength&&
				j1>=0&&j1<fftLength/2+1){

				if(!reassignObj->resultType){ // complex
					mRealArr4[i1*(fftLength/2+1)+j1]+=v1;
					mImageArr4[i1*(fftLength/2+1)+j1]+=v2;
				}
				else{ // mRealArr1 -> amp
					mRealArr4[i1*(fftLength/2+1)+j1]+=sqrtf(v1*v1+v2*v2);
				}
			}
		}
	}

	if(mRealArr5){
		memcpy(mRealArr5, mRealArr1, sizeof(float )*totalLength);
	}

	if(mImageArr5){
		memcpy(mImageArr5, mImageArr1, sizeof(float )*totalLength);
	}

}

// stft创建之后调用
static void _reassignObj_initWindowData(ReassignObj reassignObj){
	float *winArr=NULL; // S_h
	float *winDerivativeArr=NULL;
	float *winWeightArr=NULL;

	int fftLength=0;

	float *_winArr=NULL;
	float *_derArr=NULL;

	fftLength=reassignObj->fftLength;

	// 1. winArr
	winArr=__vnew(fftLength, NULL);
	_winArr=stftObj_getWindowDataArr(reassignObj->stftObj);
	memcpy(winArr, _winArr, sizeof(float )*fftLength);

	// 2. winDerivativeArr =dw(n)/dn
	winDerivativeArr=__vnew(fftLength+2, NULL);
	_derArr=__vnew(fftLength+2, NULL);
	memcpy(_derArr+1, _winArr, sizeof(float )*fftLength);
	// 'wrap' ??? 'zero'
	_derArr[0]=_derArr[fftLength];
	_derArr[fftLength+1]=_derArr[1];

	__vgradient(_derArr, fftLength+2, 1, winDerivativeArr);

	// 3. winWeightArr =n.*w(n)
	__varange(-fftLength/2, fftLength/2+1, 1, &winWeightArr);
	__vmul(winWeightArr, winArr, fftLength, winWeightArr);

	reassignObj->winArr=winArr;
	reassignObj->winDerivativeArr=winDerivativeArr; // ???
	reassignObj->winWeightArr=winWeightArr;

	free(_derArr);
}

// timeLength correlation相关数据 xcorr =>tcorr
static void _reassignObj_dealTCorrData(ReassignObj reassignObj,int dataLength){
	int timeLength=0;
	int fftLength=0;

	float *timeArr=NULL; // timeLength
	
	float *mReFreArr=NULL; // timeLength*(fftLength/2+1)
	float *mReTimeArr=NULL;

	float *mRealArr1=NULL; // S_h
	float *mImageArr1=NULL;

	float *mRealArr2=NULL; // S_dh
	float *mImageArr2=NULL;

	float *mRealArr3=NULL; // S_th
	float *mImageArr3=NULL;

	STFTObj stftObj=NULL;

	int *mTimeIndexArr=NULL;
	int *mFreIndexArr=NULL;

	int *mTempIndexArr=NULL;

	stftObj=reassignObj->stftObj;

	fftLength=reassignObj->fftLength;

	timeArr=reassignObj->timeArr;

	mReFreArr=reassignObj->mReFreArr;
	mReTimeArr=reassignObj->mReTimeArr;

	mRealArr1=reassignObj->mRealArr1;
	mImageArr1=reassignObj->mImageArr1;

	mRealArr2=reassignObj->mRealArr2;
	mImageArr2=reassignObj->mImageArr2;

	mRealArr3=reassignObj->mRealArr3;
	mImageArr3=reassignObj->mImageArr3;

	mTimeIndexArr=reassignObj->mTimeIndexArr;
	mFreIndexArr=reassignObj->mFreIndexArr;

	mTempIndexArr=reassignObj->mTempIndexArr;

	timeLength=stftObj_calTimeLength(stftObj,dataLength);
	if(reassignObj->timeLength<timeLength||
		reassignObj->timeLength>timeLength*2){
		free(mRealArr1);
		free(mImageArr1);

		free(mRealArr2);
		free(mImageArr2);

		free(mRealArr3);
		free(mImageArr3);

		free(timeArr);
		
		free(mReTimeArr);
		free(mReFreArr);

		free(mTimeIndexArr);
		free(mFreIndexArr);

		free(mTempIndexArr);

		mRealArr1=__vnew(timeLength*fftLength, NULL);
		mImageArr1=__vnew(timeLength*fftLength, NULL);

		mRealArr2=__vnew(timeLength*fftLength, NULL);
		mImageArr2=__vnew(timeLength*fftLength, NULL);

		mRealArr3=__vnew(timeLength*fftLength, NULL);
		mImageArr3=__vnew(timeLength*fftLength, NULL);

		mReTimeArr=__vnew(timeLength*(fftLength/2+1), NULL);
		mReFreArr=__vnew(timeLength*(fftLength/2+1), NULL);

		mTimeIndexArr=__vnewi(timeLength*(fftLength/2+1), NULL);
		mFreIndexArr=__vnewi(timeLength*(fftLength/2+1), NULL);

		mTempIndexArr=__vnewi(timeLength*(fftLength/2+1), NULL);

		__varange(0, timeLength, 1, &timeArr);
		for(int i=0;i<timeLength;i++){
			float _t=0;

			// isPadding True 和fftLength无关 False +fftLength/2 offset ???
			// _t=timeArr[i]*reassignObj->slideLength+(reassignObj->isPad?0:fftLength/2);
			_t=timeArr[i]*reassignObj->slideLength;
			timeArr[i]=_t/reassignObj->samplate;
		}
	}

	reassignObj->timeLength=timeLength;

	reassignObj->timeArr=timeArr;

	reassignObj->mReTimeArr=mReTimeArr;
	reassignObj->mReFreArr=mReFreArr;

	reassignObj->mRealArr1=mRealArr1;
	reassignObj->mImageArr1=mImageArr1;

	reassignObj->mRealArr2=mRealArr2;
	reassignObj->mImageArr2=mImageArr2;

	reassignObj->mRealArr3=mRealArr3;
	reassignObj->mImageArr3=mImageArr3;

	reassignObj->mTimeIndexArr=mTimeIndexArr;
	reassignObj->mFreIndexArr=mFreIndexArr;

	reassignObj->mTempIndexArr=mTempIndexArr;
}

static void _reassignObj_stft(ReassignObj reassignObj,float *dataArr,int dataLength){
	STFTObj stftObj=NULL;
	float *winArr=NULL; 

	int fftLength=0; 
	int timeLength=0;

	float *mRealArr1=NULL; // S_h
	float *mImageArr1=NULL;

	stftObj=reassignObj->stftObj;
	winArr=reassignObj->winArr;

	mRealArr1=reassignObj->mRealArr1;
	mImageArr1=reassignObj->mImageArr1;

	fftLength=reassignObj->fftLength;
	timeLength=reassignObj->timeLength;

	// reset windData & stft
	stftObj_useWindowDataArr(stftObj, winArr);
	stftObj_stft(stftObj, dataArr, dataLength, mRealArr1, mImageArr1);

	// timeLength*fftLength => timeLength*(fftLength/2+1) 
	__mccut(mRealArr1,mImageArr1,
			timeLength,fftLength,
			0,timeLength,
			0,fftLength/2+1,
			mRealArr1,mImageArr1);
}

/***
	w=w-Image(S_dh/S_h)
	t=t+Real(S_th/S_h)
****/
static void _reassignObj_reassignTimeFre(ReassignObj reassignObj,float *dataArr,int dataLength){
	ReassignType resType=Reassign_All;

	int samplate=0;
	STFTObj stftObj=NULL;

	int fftLength=0; 
	int timeLength=0;

	float *freArr=NULL;
	float *timeArr=NULL;

	float *mReFreArr=NULL;
	float *mReTimeArr=NULL;

	float *winDerivativeArr=NULL; // +1
	float *winWeightArr=NULL;

	float *mRealArr1=NULL; // S_h
	float *mImageArr1=NULL;

	float *mRealArr2=NULL; // S_dh
	float *mImageArr2=NULL;

	float *mRealArr3=NULL; // S_th
	float *mImageArr3=NULL;

	resType=reassignObj->resType;

	samplate=reassignObj->samplate;
	stftObj=reassignObj->stftObj;

	fftLength=reassignObj->fftLength;
	timeLength=reassignObj->timeLength;

	mRealArr1=reassignObj->mRealArr1;
	mImageArr1=reassignObj->mImageArr1;

	mRealArr2=reassignObj->mRealArr2;
	mImageArr2=reassignObj->mImageArr2;

	mRealArr3=reassignObj->mRealArr3;
	mImageArr3=reassignObj->mImageArr3;

	winDerivativeArr=reassignObj->winDerivativeArr;
	winWeightArr=reassignObj->winWeightArr;

	freArr=reassignObj->freArr;
	timeArr=reassignObj->timeArr;

	mReFreArr=reassignObj->mReFreArr;
	mReTimeArr=reassignObj->mReTimeArr;

	if(resType==Reassign_Fre||resType==Reassign_All){
		stftObj_useWindowDataArr(stftObj, winDerivativeArr+1); // +1
		stftObj_stft(stftObj, dataArr, dataLength, mRealArr2, mImageArr2);

		// timeLength*fftLength => timeLength*(fftLength/2+1) 
		__mccut(mRealArr2,mImageArr2,
				timeLength,fftLength,
				0,timeLength,
				0,fftLength/2+1,
				mRealArr2,mImageArr2);

	}

	if(resType==Reassign_Time||resType==Reassign_All){
		stftObj_useWindowDataArr(stftObj, winWeightArr);
		stftObj_stft(stftObj, dataArr, dataLength, mRealArr3, mImageArr3);

		// timeLength*fftLength => timeLength*(fftLength/2+1) 
		__mccut(mRealArr3,mImageArr3,
				timeLength,fftLength,
				0,timeLength,
				0,fftLength/2+1,
				mRealArr3,mImageArr3);
	}

	// mReFreArr w=w-image(S_dh/S_h)
	if(resType==Reassign_Fre||resType==Reassign_All){
		__mcdiv(mRealArr2,mImageArr2,
				mRealArr1,mImageArr1,
				timeLength,fftLength/2+1,
				mRealArr2,mImageArr2);

		__mmul_value(mImageArr2,-0.5*samplate/M_PI , timeLength, fftLength/2+1, mReFreArr);
		__madd_vector(mReFreArr, freArr,1, timeLength, fftLength/2+1, 1, NULL);
	}
	
	// mReTimeArr t=t+real(S_th/S_h)
	if(resType==Reassign_Time||resType==Reassign_All){
		__mcdiv(mRealArr3,mImageArr3,
				mRealArr1,mImageArr1,
				timeLength,fftLength/2+1,
				mRealArr3,mImageArr3);

		__mmul_value(mRealArr3, 1.0/samplate, timeLength, fftLength/2+1, mReTimeArr);
		__madd_vector(mReTimeArr, timeArr,0, timeLength, fftLength/2+1, 0, NULL);
	}

}

/***
	1 thresh & clip
	2. 默认值
****/
static void _reassignObj_filterTimeFre(ReassignObj reassignObj,int dataLength){
	int fftLength=0; 
	int timeLength=0;

	ReassignType resType;

	float *freArr=NULL;  // fftLength/2+1
	float *timeArr=NULL; // timeLength

	float *mReFreArr=NULL; // timeLength*(fftLength/2+1)
	float *mReTimeArr=NULL;

	float *mRealArr1=0; // S_h
	float *mImageArr1=0;

	int samplate=0;
	float thresh=0; // 0.001
	int mLen=0;

	float fmax=0;
	float tmax=0;

	samplate=reassignObj->samplate;
	thresh=reassignObj->thresh;
	resType=reassignObj->resType;

	fftLength=reassignObj->fftLength;
	timeLength=reassignObj->timeLength;

	freArr=reassignObj->freArr;
	timeArr=reassignObj->timeArr;

	mReFreArr=reassignObj->mReFreArr;
	mReTimeArr=reassignObj->mReTimeArr;

	mRealArr1=reassignObj->mRealArr1;
	mImageArr1=reassignObj->mImageArr1;

	mLen=fftLength/2+1;

	fmax=freArr[fftLength/2];
	tmax=timeArr[timeLength-1]; // ??? isPadding

	// debug
	// {
	// 	printf("mReFreArr is:\n");
	// 	__mdebug(mReFreArr, timeLength, fftLength/2+1, 1);
	// 	printf("\n");

	// 	printf("mReTimeArr is:\n");
	// 	__mdebug(mReTimeArr, timeLength, fftLength/2+1, 1);
	// 	printf("\n");
	// }

	// 1. thresh & clip 
	for(int i=0;i<timeLength;i++){
		for(int j=0;j<mLen;j++){
			float _power=0;
			float _r1=0;
			float _i1=0;

			float _value2=0;

			_r1=mRealArr1[i*mLen+j];
			_i1=mImageArr1[i*mLen+j];
			_power=_r1*_r1+_i1*_i1;

			if(resType==Reassign_Fre||resType==Reassign_All){ // fre
				// thresh 避免NaN
				if(_power>=thresh*thresh){
					_value2=mReFreArr[i*mLen+j];
				}
				else{
					_value2=freArr[j];
				}

				// clip
				if(_value2<0){
					_value2=0;
				}

				if(_value2>fmax){
					_value2=fmax;
				}

				mReFreArr[i*mLen+j]=_value2;
			}

			if(resType==Reassign_Time||resType==Reassign_All){ // time
				// thresh 避免NaN
				if(_power>=thresh*thresh){
					_value2=mReTimeArr[i*mLen+j];
				}
				else{
					_value2=timeArr[i];
				}

				// clip
				if(_value2<0){
					_value2=0;
				}

				if(_value2>tmax){
					_value2=tmax;
				}

				mReTimeArr[i*mLen+j]=_value2;
			}
		}
	}

	// 3. 默认值
	if(resType!=Reassign_All){
		if(resType!=Reassign_Fre){
			__mrepeat(freArr, 1, mLen, 0, timeLength, mReFreArr);
		}
		else if(resType!=Reassign_Time){
			__mrepeat(timeArr, timeLength, 1, 1, mLen, mReTimeArr);
		}
	}
}

void reassignObj_free(ReassignObj reassignObj){
	STFTObj stftObj=NULL; 

	float *freArr=NULL;  // fftLength/2+1
	float *timeArr=NULL; // timeLength

	float *mReFreArr=NULL; // timeLength*(fftLength/2+1)
	float *mReTimeArr=NULL;

	float *winArr=NULL; // S_h
	float *winDerivativeArr=NULL; // +1
	float *winWeightArr=NULL;

	// stft相关 timeLength*fftLength => timeLength*(fftLength/2+1)
	float *mRealArr1=NULL; // S_h
	float *mImageArr1=NULL;

	float *mRealArr2=NULL; // S_dh
	float *mImageArr2=NULL;

	float *mRealArr3=NULL; // S_th
	float *mImageArr3=NULL;

	int *mTimeIndexArr=NULL;
	int *mFreIndexArr=NULL;

	int *mTempIndexArr=NULL;

	if(!reassignObj){
		return;
	}

	stftObj=reassignObj->stftObj;

	freArr=reassignObj->freArr;
	timeArr=reassignObj->timeArr;

	mReFreArr=reassignObj->mReFreArr;
	mReTimeArr=reassignObj->mReTimeArr;

	winArr=reassignObj->winArr;
	winDerivativeArr=reassignObj->winDerivativeArr;
	winWeightArr=reassignObj->winWeightArr;

	mRealArr1=reassignObj->mRealArr1;
	mImageArr1=reassignObj->mImageArr1;

	mRealArr2=reassignObj->mRealArr2;
	mImageArr2=reassignObj->mImageArr2;

	mRealArr3=reassignObj->mRealArr3;
	mImageArr3=reassignObj->mImageArr3;

	mTimeIndexArr=reassignObj->mTimeIndexArr;
	mFreIndexArr=reassignObj->mFreIndexArr;

	mTempIndexArr=reassignObj->mTempIndexArr;

	stftObj_free(stftObj);

	free(freArr);
	free(timeArr);

	free(mReFreArr);
	free(mReTimeArr);

	free(winArr);
	free(winDerivativeArr);
	free(winWeightArr);

	free(mRealArr1);
	free(mImageArr1);

	free(mRealArr2);
	free(mImageArr2);

	free(mRealArr3);
	free(mImageArr3);

	free(mTimeIndexArr);
	free(mFreIndexArr);

	free(mTempIndexArr);

	free(reassignObj);
}




