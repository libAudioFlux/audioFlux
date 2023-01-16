// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/fft_algorithm.h"
#include "../dsp/dct_algorithm.h"

#include "xxcc_algorithm.h"

struct OpaqueXXCC{
	int num;
	int timeLength;

	// dct
	FFTObj fftObj; // dct accelerate
	DCTObj dctObj; // direct matrix

	float *mRealArr; // timeLength*num
	float *mImageArr;

};

static int _calRadix2(int length);

int xxccObj_new(XXCCObj *xxccObj,int num){
	int status=0;

	if(num<2){
		printf("num is error!!!\n");
		return -1;
	}

	FFTObj fftObj=NULL;
	DCTObj dctObj=NULL;
	int r=0;

	XXCCObj xc=NULL;

	xc=*xxccObj=(XXCCObj )calloc(1, sizeof(struct OpaqueXXCC ));

	r=_calRadix2(num);
	if(r){ // dct accelerate
		fftObj_new(&fftObj, r);
	}
	else{ // direct matrix
		dctObj_new(&dctObj, num,NULL);
	}

	xc->num=num;

	xc->fftObj=fftObj;
	xc->dctObj=dctObj;

	return status;
}

void xxccObj_setTimeLength(XXCCObj xxccObj,int timeLength){
	int num=0;

	float *mRealArr=NULL;
	float *mImageArr=NULL;

	num=xxccObj->num;

	mRealArr=xxccObj->mRealArr; // timeLength*num
	mImageArr=xxccObj->mImageArr;

	if(xxccObj->timeLength<timeLength||
		xxccObj->timeLength>timeLength*2){
		
		free(mRealArr);
		free(mImageArr);

		mRealArr=__vnew(timeLength*num, NULL);
		mImageArr=__vnew(timeLength*num, NULL);
	}

	xxccObj->timeLength=timeLength;

	xxccObj->mRealArr=mRealArr;
	xxccObj->mImageArr=mImageArr;
}

/***
	1. log&DCT
	2. delta/deltaDelta
****/
void xxccObj_xxcc(XXCCObj xxccObj,float *mDataArr1,int mLength,CepstralRectifyType *rectifyType,float *mDataArr2){
	FFTObj fftObj=NULL;
	DCTObj dctObj=NULL;

	CepstralRectifyType recType=CepstralRectify_Log;

	int num=0;
	int timeLength=0;

	float *mRealArr=NULL; 
	float *mImageArr=NULL;

	num=xxccObj->num;
	timeLength=xxccObj->timeLength;

	mRealArr=xxccObj->mRealArr;
	mImageArr=xxccObj->mImageArr;

	fftObj=xxccObj->fftObj;
	dctObj=xxccObj->dctObj;

	if(mLength>num){
		return;
	}

	if(rectifyType){
		recType=*rectifyType;
	}

	if(recType==CepstralRectify_CubicRoot){
		for(int i=0;i<timeLength*num;i++){
			mRealArr[i]=powf(mDataArr1[i], 1.0/3);
		}
	}
	else { // log
		for(int i=0;i<timeLength*num;i++){
			float _value=0;

			_value=mDataArr1[i];
			if(_value<1e-8){
				_value=1e-8;
			}

			mRealArr[i]=log10f(_value); // matlab canonic
		}
	}

	for(int i=0;i<timeLength;i++){
		if(fftObj){
			fftObj_dct(fftObj, mRealArr+i*num, mImageArr+i*num, 1);
		}
		else{
			dctObj_dct(dctObj, mRealArr+i*num,1,mImageArr+i*num);
		}
	}

	for(int i=0;i<timeLength;i++){
		for(int j=0;j<mLength;j++){
			mDataArr2[i*mLength+j]=mImageArr[i*num+j];
		}
	}
}

/***
	mfcc standard/xxcc standard
	mDataArr1 'mag' mag||power
	ccNum<=num
	[ccNum,delta,deltaDelta] vector
	deltaWindowLength 9(defalut); must odd>=3
	energyType Repalce; Append ccNum+1
	rectifyType Log
	return timeLength*[ccNum,delta,deltaDelta]
****/
void xxccObj_xxccStandard(XXCCObj xxccObj,float *mDataArr1,int mLength,float *energyArr,
						int *deltaWindowLength,CepstralEnergyType *energyType,CepstralRectifyType *rectifyType,
						float *mCoeArr,float *mDeltaArr1,float *mDeltaArr2){
	int dLen=9;
	CepstralEnergyType eType=CepstralEnergy_Replace;
	CepstralRectifyType rType=CepstralRectify_Log;

	FFTObj fftObj=NULL;
	DCTObj dctObj=NULL;

	int num=0;
	int timeLength=0;

	float *mRealArr=NULL; 
	float *mImageArr=NULL;

	float energy=0;

	float *vArr1=NULL;
	float *vArr2=NULL;

	num=xxccObj->num;
	timeLength=xxccObj->timeLength;

	mRealArr=xxccObj->mRealArr;
	mImageArr=xxccObj->mImageArr;

	fftObj=xxccObj->fftObj;
	dctObj=xxccObj->dctObj;

	if(mLength>num){
		return;
	}

	if(deltaWindowLength){
		if(*deltaWindowLength>=3&&(*deltaWindowLength)%2==1){
			dLen=*deltaWindowLength;
		}
	}

	if(energyType){
		eType=*energyType;
	}

	if(rectifyType){
		rType=*rectifyType;
	}

	if(rType==CepstralRectify_CubicRoot){
		for(int i=0;i<timeLength*num;i++){
			mRealArr[i]=powf(mDataArr1[i], 1.0/3);
		}
	}
	else { // log
		for(int i=0;i<timeLength*num;i++){
			float _value=0;

			_value=mDataArr1[i];
			if(_value<1e-8){
				_value=1e-8;
			}

			mRealArr[i]=log10f(_value); // matlab canonic
		}
	}

	for(int i=0;i<timeLength;i++){
		if(fftObj){
			fftObj_dct(fftObj, mRealArr+i*num, mImageArr+i*num, 1);
		}
		else{
			dctObj_dct(dctObj, mRealArr+i*num,1,mImageArr+i*num);
		}
	}

	// 1. coef
	for(int i=0;i<timeLength;i++){
		if(eType!=CepstralEnergy_Ignore){
			if(energyArr[i]<1e-8){
				energy=logf(1e-8);
			}
			else{
				energy=logf(energyArr[i]);
			}
		}

		for(int j=0;j<mLength;j++){
			if(eType==CepstralEnergy_Replace){
				if(!j){
					mCoeArr[i*mLength+j]=energy;
				}
				else{
					mCoeArr[i*mLength+j]=mImageArr[i*num+j];
				}
			}
			else if(eType==CepstralEnergy_Append){
				if(!j){
					mCoeArr[i*(mLength+1)+j]=energy;
				}
				
				mCoeArr[i*(mLength+1)+(j+1)]=mImageArr[i*num+j];
			}
			else{
				mCoeArr[i*mLength+j]=mImageArr[i*num+j];
			}
		}
	}

	// 2. delta&deltaDelta
	if(eType==CepstralEnergy_Append){
		mLength+=1;
	}

	vArr1=__vnew(mLength, NULL);
	vArr2=__vnew(mLength, NULL);
	for(int i=0;i<timeLength;i++){
		util_delta(mCoeArr+i*mLength,mLength,dLen,vArr1);
		memcpy(mDeltaArr1+i*mLength,vArr1,sizeof(float )*mLength);

		util_delta(vArr1,mLength,dLen,vArr2);
		memcpy(mDeltaArr2+i*mLength,vArr2,sizeof(float )*mLength);

		memset(vArr1, 0, sizeof(float )*mLength);
		memset(vArr2, 0, sizeof(float )*mLength);
	}

	free(vArr1);
	free(vArr2);
}

void xxccObj_free(XXCCObj xxccObj){
	FFTObj fftObj=NULL; // dct accelerate
	DCTObj dctObj=NULL; // direct matrix

	float *mRealArr=NULL; // timeLength*num
	float *mImageArr=NULL;

	if(xxccObj){
		fftObj=xxccObj->fftObj;
		dctObj=xxccObj->dctObj;

		mRealArr=xxccObj->mRealArr;
		mImageArr=xxccObj->mImageArr;

		fftObj_free(fftObj);
		dctObj_free(dctObj);

		free(mRealArr);
		free(mImageArr);

		free(xxccObj);
	}
}

static int _calRadix2(int length){
	int r=0;
	int m=0;

	while(1){
		m=length%2;
		if(!m){ //
			r++;
			length=length/2;
			if(length==1){
				break;
			}
		}
		else{
			r=0;
			break;
		}
	}

	return r;
}













