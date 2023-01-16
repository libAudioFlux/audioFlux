// clang 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "synsq_algorithm.h"

struct OpaqueSynsq{
	int samplate;

	int num; // fre
	int fftLength; // time

	int order; // >=1

	int *mFreIndexArr; // fre*time
	int *mTimeIndexArr;

	int *mTempIndexArr;

	float *mArr1;
	float *mArr2;

	float *vArr1; // num

	float thresh;
};

// for /bark/erb
static int __arr_roundIndex(float *arr,int length,float value);

int synsqObj_new(SynsqObj *synsqObj,int num,int radix2Exp,
				int *samplate,int *order,
				float *thresh){
	int status=0;
	SynsqObj ss=NULL;

	int fftLength=0;
	int _samplate=32000;
	int _order=1;

	int *mFreIndexArr=NULL; // fre*time
	int *mTimeIndexArr=NULL;

	int *mTempIndexArr=NULL;

	float *mArr1=NULL;
	float *mArr2=NULL;

	float *vArr1=NULL;

	float _thresh=0.001;

	ss=*synsqObj=(SynsqObj )calloc(1, sizeof(struct OpaqueSynsq ));

	fftLength=1<<radix2Exp;
	if(samplate){
		if(*samplate>0&&*samplate<196000){
			_samplate=*samplate;
		}
	}

	if(order){
		if(*order>1){
			_order=*order;
		}
	}

	if(thresh){
		if(*thresh>1){
			_thresh=*thresh;
		}
	}

	mFreIndexArr=__vnewi(num*fftLength, NULL);
	mTimeIndexArr=__vnewi(num*fftLength, NULL);

	mTempIndexArr=__vnewi(num*fftLength, NULL);

	mArr1=__vnew(num*fftLength, NULL);
	mArr2=__vnew(num*fftLength, NULL);

	vArr1=__vnew(num, NULL);

	mTimeIndexArr=__vnewi(num*fftLength, NULL);
	for(int i=0;i<fftLength;i++){
		mTimeIndexArr[i]=i;
	}
	for(int i=1;i<num;i++){
		memcpy(mTimeIndexArr+i*fftLength, mTimeIndexArr, sizeof(int )*fftLength);
	}

	ss->samplate=_samplate;

	ss->num=num;
	ss->fftLength=fftLength;

	ss->order=_order;

	ss->mFreIndexArr=mFreIndexArr;
	ss->mTimeIndexArr=mTimeIndexArr;

	ss->mTempIndexArr=mTempIndexArr;

	ss->mArr1=mArr1;
	ss->mArr2=mArr2;

	ss->vArr1=vArr1;

	ss->thresh=_thresh;

	return status;
}

/***
	1. angle
	2. unwrap
	3. diff
	4. index
	5. order
	6. result
****/
void synsqObj_synsq(SynsqObj synsqObj,float *freArr,
					SpectralFilterBankScaleType scaleType,
					float *mRealArr1,float *mImageArr1,
					float *mRealArr2,float *mImageArr2){

	int samplate=0;

	int num=0; // fre
	int fftLength=0; // time

	int order=1; // >=1

	int *mFreIndexArr=NULL; // fre*time
	int *mTimeIndexArr=NULL;

	int *mTempIndexArr=NULL;

	float *mArr1=NULL;
	float *mArr2=NULL;	

	float *vArr1=NULL;

	float thresh=0;
	
	float fmin=0;
	float fmax=0;

	samplate=synsqObj->samplate;

	num=synsqObj->num;
	fftLength=synsqObj->fftLength;

	order=synsqObj->order;

	mFreIndexArr=synsqObj->mFreIndexArr;
	mTimeIndexArr=synsqObj->mTimeIndexArr;

	mTempIndexArr=synsqObj->mTempIndexArr;

	mArr1=synsqObj->mArr1;
	mArr2=synsqObj->mArr2;

	vArr1=synsqObj->vArr1;

	thresh=synsqObj->thresh;

	if(scaleType>SpectralFilterBankScale_Log){
		printf("scaleType is error!\n");
		return ;
	}

	// 1. angle
	for(int i=0;i<num*fftLength;i++){
		mArr1[i]=atan2f(mRealArr1[i], mImageArr1[i]);
	}

	// 2. unwarp
	__munwrap(mArr1, num,fftLength, 1);

	// 3. diff
	__mdiff2(mArr1, num, fftLength, 1, NULL, mArr2);
	for(int i=0;i<num;i++){
		mArr2[i*fftLength+fftLength-1]=mArr2[i*fftLength+fftLength-2];
	}
	__vdiv_value(mArr2, 2*M_PI, num*fftLength, NULL);
	
	// 4. index
	if(scaleType==SpectralFilterBankScale_Octave||
		scaleType==SpectralFilterBankScale_Log){ // log 

		fmin=freArr[0]/samplate;
		fmax=freArr[num-1]/samplate;
		for(int i=0;i<num*fftLength;i++){ // floorf ???
			mFreIndexArr[i]=roundf((log2f(fabsf(mArr2[i]))-log2f(fmin))*num/(log2f(fmax)-log2f(fmin)));
		}
	}
	else if(scaleType==SpectralFilterBankScale_Linear||
		scaleType==SpectralFilterBankScale_Linspace){ // linear 
		
		fmin=freArr[0]/samplate;
		fmax=freArr[num-1]/samplate;
		for(int i=0;i<num*fftLength;i++){ // floorf ???
			mFreIndexArr[i]=roundf(fabsf((mArr2[i])-fmin)*num/(fmax-fmin));
		}
	}
	else{ // mel/bark/erb ???
		__vdiv_value(freArr, samplate, num, vArr1);
		for(int i=0;i<num*fftLength;i++){ // floorf ???
			mFreIndexArr[i]=__arr_roundIndex(vArr1, num, fabsf(mArr2[i]));
		}
	}

	// 5. order
	if(order>1){
		int v=0;

		memset(mTempIndexArr, 0, sizeof(int )*num*fftLength);
		for(int k=0;k<order-1;k++){

			for(int i=0;i<fftLength;i++){
				for(int j=0;j<num;j++){
					
					v=mFreIndexArr[i*num+j];
					if(v>=0&&v<num){
						mTempIndexArr[i*num+j]=mFreIndexArr[i*num+v];
					}
				}
			}

			memcpy(mFreIndexArr, mTempIndexArr, sizeof(int )*num*fftLength);
		}
	}

	// 6. result
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

				mRealArr2[i1*fftLength+j1]+=v1;
				mImageArr2[i1*fftLength+j1]+=v2;
			}
		}
	}

}

void synsqObj_free(SynsqObj synsqObj){
	int *mFreIndexArr=NULL; // fre*time
	int *mTimeIndexArr=NULL;

	int *mTempIndexArr=NULL;

	float *mArr1=NULL;
	float *mArr2=NULL;

	float *vArr1=NULL;

	if(synsqObj){
		mFreIndexArr=synsqObj->mFreIndexArr;
		mTimeIndexArr=synsqObj->mTimeIndexArr;

		mTempIndexArr=synsqObj->mTempIndexArr;

		mArr1=synsqObj->mArr1;
		mArr2=synsqObj->mArr2;

		vArr1=synsqObj->vArr1;

		free(mFreIndexArr);
		free(mTimeIndexArr);

		free(mTempIndexArr);

		free(mArr1);
		free(mArr2);

		free(vArr1);

		free(synsqObj);
	}
}

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









