// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "cwt_algorithm.h"
#include "wsst_algorithm.h"

struct OpaqueWSST{
	CWTObj cwtObj;

	int num;
	int fftLength; // 1<<radix2Exp
	
	// fre*time=num*fftLength
	float *mRealArr1;  // cwt
	float *mImageArr1;

	float *mRealArr2; //cwtDet
	float *mImageArr2;

	float *mRealArr3; // phase ->mImageArr3
	float *mImageArr3;

	// rearrage相关 num*fftLength
	int *mFreIndexArr;
	int *mTimeIndexArr;

	int *mTempIndexArr;

	float *vArr1; // num

	float thresh; // 0.001
	int samplate;

	WaveletContinueType waveletType;
	SpectralFilterBankScaleType scaleType;

	int order; // >=1 
};

// for /bark/erb
static int __arr_roundIndex(float *arr,int length,float value);

/***
	thresh >=0 default 0.0001
	waveletType 'morlet'
	scaleType 'log'

	'morse' gamma 3 beta 20
	'morlet' gamma 6 beta 2
	'bump' gamma 5 beta 0.6

	'paul' gamma 4
	'dog' gamma 2 beta 2; must even
	'mexican' beta 2 
****/
int wsstObj_new(WSSTObj *wsstObj,int num,int radix2Exp,
			 	int *samplate,float *lowFre,float *highFre,int *binPerOctave,
			 	WaveletContinueType *waveletType,
			 	SpectralFilterBankScaleType *scaleType,
			 	float *gamma,float *beta,
			 	float *thresh,
			 	int *isPadding){
	int status=0;
	WSSTObj wsst=NULL;

	int fftLength=0;
	CWTObj cwtObj=NULL;

	int *mTimeIndexArr=NULL;

	float _thresh=0.001;
	int _samplate=32000;

	WaveletContinueType _waveletType=WaveletContinue_Morlet;
	SpectralFilterBankScaleType _scaleType=SpectralFilterBankScale_Octave;

	if(thresh){
		if(*thresh>=0){
			_thresh=*thresh;
		}
	}

	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
		}
	}

	if(waveletType){
		_waveletType=*waveletType;
	}

	if(scaleType){
		_scaleType=*scaleType;
		if(_scaleType>SpectralFilterBankScale_Log){
			printf("scaleType is error!\n");
			return 1;
		}
	}

	fftLength=1<<radix2Exp;
	wsst=*wsstObj=(WSSTObj )calloc(1, sizeof(struct OpaqueWSST ));

	status=cwtObj_new(&cwtObj,num,radix2Exp,
			samplate,lowFre,highFre,binPerOctave,
			&_waveletType,
			&_scaleType,
			gamma,beta,
			isPadding);
	cwtObj_enableDet(cwtObj, 1);

	mTimeIndexArr=__vnewi(num*fftLength, NULL);
	for(int i=0;i<fftLength;i++){
		mTimeIndexArr[i]=i;
	}
	for(int i=1;i<num;i++){
		memcpy(mTimeIndexArr+i*fftLength, mTimeIndexArr, sizeof(int )*fftLength);
	}

	wsst->cwtObj=cwtObj;

	wsst->num=num;
	wsst->fftLength=fftLength;

	wsst->mRealArr1=__vnew(num*fftLength, NULL);
	wsst->mImageArr1=__vnew(num*fftLength, NULL);

	wsst->mRealArr2=__vnew(num*fftLength, NULL);
	wsst->mImageArr2=__vnew(num*fftLength, NULL);

	wsst->mRealArr3=__vnew(num*fftLength, NULL);
	wsst->mImageArr3=__vnew(num*fftLength, NULL);

	wsst->mFreIndexArr=__vnewi(num*fftLength, NULL);
	wsst->mTimeIndexArr=mTimeIndexArr;

	wsst->vArr1=__vnew(num, NULL);;

	wsst->thresh=_thresh;
	wsst->samplate=_samplate;

	wsst->waveletType=_waveletType;
	wsst->scaleType=_scaleType;

	return status;
}

float *wsstObj_getFreBandArr(WSSTObj wsstObj){

	return cwtObj_getFreBandArr(wsstObj->cwtObj);
}

int *wsstObj_getBinBandArr(WSSTObj wsstObj){

	return cwtObj_getBinBandArr(wsstObj->cwtObj);
}

// order >=1
void wsstObj_setOrder(WSSTObj wsstObj,int order){

	wsstObj->order=order;
}

// mRealArr5&mImageArr5 may NULL
/***
	1. phase
	2. rearrage
****/
void wsstObj_wsst(WSSTObj wsstObj,float *dataArr,
				float *mRealArr4,float *mImageArr4,
				float *mRealArr5,float *mImageArr5){
	CWTObj cwtObj=NULL;

	int num=0;
	int fftLength=0; // 1<<radix2Exp
	
	// fre*time=num*fftLength
	float *mRealArr1=NULL;  // cwt
	float *mImageArr1=NULL;

	float *mRealArr2=NULL; //cwtDet
	float *mImageArr2=NULL;

	float *mRealArr3=NULL; // phase
	float *mImageArr3=NULL;

	int *mFreIndexArr=NULL; 
	int *mTimeIndexArr=NULL;

	int *mTempIndexArr=NULL;

	float *vArr1=NULL;

	float thresh=0;
	int samplate=0;

	int totalLength=0;

	float fmin=0;
	float fmax=0;

	float *freArr=NULL;
	
	int order=1;

	cwtObj=wsstObj->cwtObj;

	num=wsstObj->num;
	fftLength=wsstObj->fftLength;

	mRealArr1=wsstObj->mRealArr1;
	mImageArr1=wsstObj->mImageArr1;

	mRealArr2=wsstObj->mRealArr2;
	mImageArr2=wsstObj->mImageArr2;

	mRealArr3=wsstObj->mRealArr3;
	mImageArr3=wsstObj->mImageArr3;

	mFreIndexArr=wsstObj->mFreIndexArr;
	mTimeIndexArr=wsstObj->mTimeIndexArr;

	mTempIndexArr=wsstObj->mTempIndexArr;
	
	vArr1=wsstObj->vArr1;

	thresh=wsstObj->thresh;
	samplate=wsstObj->samplate;

	totalLength=num*fftLength;

	order=wsstObj->order;

	freArr=cwtObj_getFreBandArr(cwtObj);

	// 1. phase
	cwtObj_cwt(cwtObj, dataArr, mRealArr1, mImageArr1);
	cwtObj_cwtDet(cwtObj, NULL, mRealArr2, mImageArr2);

	__mcdiv(mRealArr2, mImageArr2, mRealArr1, mImageArr1, num, fftLength, mRealArr3, mImageArr3);
	__vdiv_value(mImageArr3, 2*M_PI, num*fftLength, NULL);

	// debug
	{
		// printf("cwt is :\n");
		// __mcdebug(mRealArr1, mImageArr1, num, fftLength, 1);
		// printf("\n");

		// printf("dcwt is :\n");
		// __mcdebug(mRealArr2, mImageArr2, num, fftLength, 1);
		// printf("\n");

		// printf("phase is :\n");
		// __mdebug(mImageArr3, num, fftLength, 1);
		// printf("\n");

		// __vdebug(freArr, num, 0);
	}

	// 2. rearrage
	if(wsstObj->scaleType==SpectralFilterBankScale_Octave||
		wsstObj->scaleType==SpectralFilterBankScale_Log){ // log 
		
		fmin=freArr[0]/samplate;
		fmax=freArr[num-1]/samplate;
		for(int i=0;i<totalLength;i++){ // floorf ???
			mFreIndexArr[i]=roundf((log2f(fabsf(mImageArr3[i]))-log2f(fmin))*num/(log2f(fmax)-log2f(fmin)));
		}
	}
	else if(wsstObj->scaleType==SpectralFilterBankScale_Linear||
		wsstObj->scaleType==SpectralFilterBankScale_Linspace){ // linear 
		
		fmin=freArr[0]/samplate;
		fmax=freArr[num-1]/samplate;
		for(int i=0;i<totalLength;i++){ // floorf ???
			mFreIndexArr[i]=roundf(fabsf((mImageArr3[i])-fmin)*num/(fmax-fmin));
		}
	}
	else{ // mel/bark/erb ???
		__vdiv_value(freArr, samplate, num, vArr1);
		for(int i=0;i<totalLength;i++){ // floorf ???
			mFreIndexArr[i]=__arr_roundIndex(vArr1, num, fabsf(mImageArr3[i]));
		}
	}

	// order
	if(order>1){
		int v=0;

		memset(mTempIndexArr,0,sizeof(int )*num*fftLength);
		for(int k=0;k<order-1;k++){

			for(int i=0;i<fftLength;i++){
				for(int j=0;j<num;j++){
					
					v=mFreIndexArr[i*num+j];
					if(v>=0&&v<num){
						mTempIndexArr[i*num+j]=mFreIndexArr[i*num+v];
					}
				}
			}

			memcpy(mFreIndexArr, mTempIndexArr, sizeof(int )*fftLength*num);
			// memset(mTempIndexArr,0,sizeof(int )*fftLength*num);
		}
	}

	for(int i=0;i<num;i++){
		int i1=0;
		int j1=0;

		float v1=0;
		float v2=0;

		for(int j=0;j<fftLength;j++){
			i1=mFreIndexArr[i*fftLength+j];
			j1=mTimeIndexArr[i*fftLength+j];

			v1=mRealArr1[i*fftLength+j];
			v2=mImageArr1[i*fftLength+j];

			if(i1>=0&&i1<num&&
				v1*v1+v2*v2>thresh*thresh){

				mRealArr4[i1*fftLength+j1]+=v1;
				mImageArr4[i1*fftLength+j1]+=v2;
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

// for /bark/erb
/***
	arr asc
****/
static int __arr_roundIndex(float *arr,int length,float value){
	int index=-1;

	float left=0;
	float right=0;

	for(int i=0;i<length-1;i++){
		if(fabsf(value)>=arr[i]&&
			fabsf(value)<arr[i+1]){

			left=fabsf(value)-arr[i];
			right=arr[i+1]-fabsf(value);

			if(left<right){
				index=i;
			}
			else{
				index=i+1;
			}
			
			break;
		}
	}

	return index;
}

void wsstObj_free(WSSTObj wsstObj){
	CWTObj cwtObj=NULL;

	float *mRealArr1=NULL;  // cwt
	float *mImageArr1=NULL;

	float *mRealArr2=NULL; //cwtDet
	float *mImageArr2=NULL;

	float *mRealArr3=NULL; // phase
	float *mImageArr3=NULL;

	int *mFreIndexArr=NULL;
	int *mTimeIndexArr=NULL;

	int *mTempIndexArr=NULL;

	float *vArr1=NULL; // num

	if(wsstObj){
		cwtObj=wsstObj->cwtObj;

		mRealArr1=wsstObj->mRealArr1;
		mImageArr1=wsstObj->mImageArr1;

		mRealArr2=wsstObj->mRealArr2;
		mImageArr2=wsstObj->mImageArr2;

		mRealArr3=wsstObj->mRealArr3;
		mImageArr3=wsstObj->mImageArr3;

		mFreIndexArr=wsstObj->mFreIndexArr;
		mTimeIndexArr=wsstObj->mTimeIndexArr;

		mTempIndexArr=wsstObj->mTempIndexArr;

		vArr1=wsstObj->vArr1;

		cwtObj_free(cwtObj);

		free(mRealArr1);
		free(mImageArr1);

		free(mRealArr2);
		free(mImageArr2);

		free(mRealArr3);
		free(mImageArr3);

		free(mFreIndexArr);
		free(mTimeIndexArr);

		free(mTempIndexArr);
		free(vArr1);

		free(wsstObj);
	}
}







