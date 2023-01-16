// 

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "viterbi.h"
#include "hmm.h"

struct OpaqueHMM{
	// λ=(π,A,B)
	float *piArr; // sLength
	float *mTransitionArr; // sLength*nLength
	float *mEmmissionArr; // sLength*nLength

	int sLength;
	int nLength;

	// forward/backward -->train cache
	float *mAlphaArr; // tLength*sLength
	float *mBetaArr; // tLength*sLength

	int tLength;

	int isDebug;
};

static void __forward(float *piArr,float *mTransitionArr,float *mEmmissionArr,
					int sLength,int nLength,
					float *mAlphaArr,int *oArr,
					int tLength);
static void __backward(float *piArr,float *mTransitionArr,float *mEmmissionArr,
					int sLength,int nLength,
					float *mBetaArr,int *oArr,
					int tLength);

static int __distribute(float *arr,int length);

static float __hmmObj_calGamma(HMMObj hmmObj,int i,int t);
static float __hmmObj_calKsi(HMMObj hmmObj,int *oArr,int i,int j,int t);

static void __hmmObj_init(HMMObj hmmObj);

static float __hmmObj_predict(HMMObj hmmObj,int *oArr,int tLength,int type);

int hmmObj_new(HMMObj *hmmObj,int sLength,int nLength){
	int status=0;

	float *piArr=NULL; // λ=(π,A,B)
	float *mTransitionArr=NULL;
	float *mEmmissionArr=NULL;

	HMMObj hmm=NULL;

	hmm=*hmmObj=(HMMObj )calloc(1, sizeof(struct OpaqueHMM ));

	piArr=__vnew(sLength, NULL);
	mTransitionArr=__vnew(sLength*sLength, NULL);
	mEmmissionArr=__vnew(sLength*nLength, NULL);

	hmm->piArr=piArr;
	hmm->mTransitionArr=mTransitionArr;
	hmm->mEmmissionArr=mEmmissionArr;

	hmm->sLength=sLength;
	hmm->nLength=nLength;

	__hmmObj_init(hmm);

	return status;
}

void hmmObj_init(HMMObj hmmObj,float *piArr,float *mAArr,float *mBArr){
	float sum=0;
	float *arr1=NULL;

	arr1=__vnew(hmmObj->sLength, NULL);

	if(piArr){
		sum=__vsum(piArr, hmmObj->sLength);
		if(sum==1.0){
			free(hmmObj->piArr);
			hmmObj->piArr=piArr;
		}
		else{
			printf("init pi is fail !!!\n");
		}
	}

	if(mAArr){
		__msum(mAArr, hmmObj->sLength, hmmObj->sLength, 1, arr1);
		sum=__vsum(arr1, hmmObj->sLength);
		if(sum==hmmObj->sLength*1.0){
			free(hmmObj->mTransitionArr);
			hmmObj->mTransitionArr=mAArr;
		}
		else{
			printf("init A matrix is fail !!!\n");
		}
	}

	if(mBArr){
		__msum(mBArr, hmmObj->sLength, hmmObj->nLength, 1, arr1);
		sum=__vsum(arr1, hmmObj->sLength);
		if(sum==hmmObj->sLength*1.0){
			free(hmmObj->mEmmissionArr);
			hmmObj->mEmmissionArr=mBArr;
		}
		else{
			printf("init B matrix is fail !!!\n");
		}
	}
}

void hmmObj_generate(HMMObj hmmObj,int tLength,int *oArr,int *sArr){
	int index1=0;
	int index2=0;

	float *piArr=NULL;
	float *mTransitionArr=NULL;
	float *mEmmissionArr=NULL;

	int sLength=0;
	int nLength=0;

	time_t t;	
	srand((unsigned) time(&t));

	piArr=hmmObj->piArr;
	mTransitionArr=hmmObj->mTransitionArr;
	mEmmissionArr=hmmObj->mEmmissionArr;

	sLength=hmmObj->sLength;
	nLength=hmmObj->nLength;

	index1=__distribute(piArr, sLength);
	sArr[0]=index1;
	if(hmmObj->isDebug){
		if(index1>=sLength){
			printf("sArr[%d] --> %d is error!!!\n",0,index1);
		}
	}
	
	for(int i=1;i<tLength;i++){
		index2=sArr[i-1];
		index1=__distribute(mTransitionArr+index2*sLength, sLength);

		sArr[i]=index1;
		if(hmmObj->isDebug){
			if(index1>=sLength){
				printf("sArr[%d] --> %d is error!!!\n",i,index1);
			}
		}
	}

	for(int i=0;i<tLength;i++){
		index2=sArr[i];
		index1=__distribute(mEmmissionArr+index2*nLength, nLength);

		oArr[i]=index1;
		if(hmmObj->isDebug){
			if(index1>=sLength){
				printf("oArr[%d] --> %d is error!!!\n",i,index1);
			}
		}
	}
}

void hmmObj_enableDebug(HMMObj hmmObj,int isDebug){

	hmmObj->isDebug=1;
}

/***
	Baum-Weich
	maxIter 100 error 1e-3
****/
void hmmObj_train(HMMObj hmmObj,int *oArr,int tLength,int *maxIter,float *error){
	float *vArr1=NULL;
	float *mAArr1=NULL;
	float *mBArr1=NULL;

	float *vArr2=NULL; // det
	float *mAArr2=NULL;
	float *mBArr2=NULL;

	int sLength=0;
	int nLength=0;

	float num=0;
	float den=0;

	int _maxIter=100;
	float _error=1e-3;

	float aError=0;
	float bError=0;
	float pError=0;

	sLength=hmmObj->sLength;
	nLength=hmmObj->nLength;

	vArr1=__vnew(sLength, NULL);
	mAArr1=__vnew(sLength*sLength, NULL);
	mBArr1=__vnew(sLength*nLength, NULL);

	vArr2=__vnew(sLength, NULL);
	mAArr2=__vnew(sLength*sLength, NULL);
	mBArr2=__vnew(sLength*nLength, NULL);

	if(maxIter){
		_maxIter=*maxIter;
	}

	if(error){
		_error=*error;
	}

	if(!hmmObj->mAlphaArr){
		hmmObj->mAlphaArr=__vnew(tLength*sLength, NULL);
	}
	else{
		if(tLength>hmmObj->tLength){
			free(hmmObj->mAlphaArr);
			hmmObj->mAlphaArr=__vnew(tLength*sLength, NULL);
		}
	}

	if(!hmmObj->mBetaArr){
		hmmObj->mBetaArr=__vnew(tLength*sLength, NULL);
	}
	else{
		if(tLength>hmmObj->tLength){
			free(hmmObj->mBetaArr);
			hmmObj->mBetaArr=__vnew(tLength*sLength, NULL);
		}
	}

	hmmObj->tLength=tLength;
	for(int i=0;i<_maxIter;i++){
		// 1. forward/backward --> alpha/beta matrix --> cal A/B hidden 
		memset(hmmObj->mAlphaArr, 0, sizeof(float )*tLength*sLength);
		memset(hmmObj->mBetaArr, 0, sizeof(float )*tLength*sLength);

		__forward(hmmObj->piArr, hmmObj->mTransitionArr, hmmObj->mEmmissionArr, 
				sLength, nLength, 
				hmmObj->mAlphaArr, oArr,
				tLength);

		__backward(hmmObj->piArr, hmmObj->mTransitionArr, hmmObj->mEmmissionArr, 
				sLength, nLength, 
				hmmObj->mBetaArr, oArr, 
				tLength);

		/***
			gamma γ的重复计算问题 ??? mGammaArr
			E-step/forward batch ??? E-step/forward loop???
		****/
		// 2. update A/B/pi --> cal Eξ/Eγ
		for(int i=0;i<sLength;i++){ // A sLength*sLength
			for(int j=0;j<sLength;j++){
				num=0;
				den=0;

				for(int k=0;k<tLength-1;k++){
					num+=__hmmObj_calKsi(hmmObj,oArr,i,j,k);
					den+=__hmmObj_calGamma(hmmObj, i, k);
				}

				mAArr1[i*sLength+j]=num/(den);
			}
		}

		for(int i=0;i<sLength;i++){ // B sLength*nLength
			for(int j=0;j<nLength;j++){
				num=0;
				den=0;

				for(int k=0;k<tLength;k++){
					float _value=0;

					_value=__hmmObj_calGamma(hmmObj, i, k);
					if(j==oArr[k]){
						num+=_value;
					}
					den+=_value;
				}

				mBArr1[i*nLength+j]=num/(den);
			}
		}

		for(int i=0;i<sLength;i++){ // pi
			vArr1[i]=__hmmObj_calGamma(hmmObj, i, 0);
		}

		// check error
		__vsub(hmmObj->piArr, vArr1, sLength, vArr2);
		__vsub(hmmObj->mTransitionArr, mAArr1, sLength*sLength, mAArr2);
		__vsub(hmmObj->mEmmissionArr, mBArr1, sLength*nLength, mBArr2);

		pError=__vnorm(vArr2, sLength);
		aError=__vnorm(mAArr2, sLength*sLength);
		bError=__vnorm(mBArr2, sLength*nLength);

		memcpy(hmmObj->mTransitionArr, mAArr1, sizeof(float )*sLength*sLength);
		memcpy(hmmObj->mEmmissionArr, mBArr1, sizeof(float )*sLength*nLength);
		memcpy(hmmObj->piArr, vArr1, sizeof(float )*sLength);

		if(hmmObj->isDebug){
			printf("iter %d --> A error: %f, B error: %f, pi error: %f \n",i,aError,bError,pError);

			__mdebug(hmmObj->mTransitionArr, sLength, sLength, 0);
			printf("\n");

			__mdebug(hmmObj->mEmmissionArr, sLength, nLength, 0);
			printf("\n");

			__vdebug(hmmObj->piArr, sLength, 0);
			printf("\n");
		}
		

		if(aError<=_error&&
			bError<=_error&&
			pError<=_error){

			break;
		}
	}

	free(vArr1);
	free(mAArr1);
	free(mBArr1);

	free(vArr2);
	free(mAArr2);
	free(mBArr2);
}

/***
	forward/backward
****/
float hmmObj_predict(HMMObj hmmObj,int *oArr,int tLength){
	int type=0; // forward
	float pValue=0;

	pValue=__hmmObj_predict(hmmObj,oArr,tLength,type);

	return pValue;
}

/***
	viterbi
****/
float hmmObj_decode(HMMObj hmmObj,int *oArr,int tLength,int *sArr,float *mProbArr){
	float *piArr=NULL; // sLength
	float *mTransitionArr=NULL; // sLength*nLength
	float *mEmmissionArr=NULL; // sLength*nLength

	int sLength=0;
	int nLength=0;

	float prob=0;

	piArr=hmmObj->piArr;
	mTransitionArr=hmmObj->mTransitionArr;
	mEmmissionArr=hmmObj->mEmmissionArr;

	sLength=hmmObj->sLength;
	nLength=hmmObj->nLength;

	prob=viterbi(piArr,mTransitionArr,mEmmissionArr,
			sLength,nLength,
			oArr,tLength,NULL,
			sArr,mProbArr,NULL);

	return prob;
}

void hmmObj_free(HMMObj hmmObj){
	float *piArr=NULL; // λ=(π,A,B)
	float *mTransitionArr=NULL;
	float *mEmmissionArr=NULL;

	float *mAlphaArr=NULL; // forward/backward
	float *mBetaArr=NULL;

	if(hmmObj){
		piArr=hmmObj->piArr;
		mTransitionArr=hmmObj->mTransitionArr;
		mEmmissionArr=hmmObj->mEmmissionArr;

		mAlphaArr=hmmObj->mAlphaArr;
		mBetaArr=hmmObj->mBetaArr;

		free(piArr);
		free(mTransitionArr);
		free(mEmmissionArr);

		free(mAlphaArr);
		free(mBetaArr);

		free(hmmObj);
	}
}

static void __hmmObj_init(HMMObj hmmObj){
	float *piArr=NULL; // λ=(π,A,B)
	float *mTransitionArr=NULL;
	float *mEmmissionArr=NULL;
	
	int sLength=0;	
	int nLength=0;

	int *rArr=NULL;
	int rSum=0;
	
	time_t t;	

	piArr=hmmObj->piArr;
	mTransitionArr=hmmObj->mTransitionArr;
	mEmmissionArr=hmmObj->mEmmissionArr;

	sLength=hmmObj->sLength;
	nLength=hmmObj->nLength;

	rArr=__vnewi(sLength, NULL);
	srand((unsigned) time(&t));

	// A 
	for(int i=0;i<sLength;i++){
		for(int k=0;k<sLength;k++){
			rArr[k]=rand()%100;
		}
		rSum=__vsumi(rArr, sLength);

		for(int j=0;j<sLength;j++){
			mTransitionArr[i*sLength+j]=1.0*rArr[j]/rSum;
		}
	}

	// pi
	for(int k=0;k<sLength;k++){
		rArr[k]=rand()%100;
	}
	rSum=__vsumi(rArr, sLength);

	for(int i=0;i<sLength;i++){
		piArr[i]=1.0*rArr[i]/rSum;
		// piArr[i]=1.0/sLength;
	}

	// B
	free(rArr);
	rArr=__vnewi(nLength, NULL);
	for(int i=0;i<sLength;i++){
		for(int k=0;k<nLength;k++){
			rArr[k]=rand()%100;
		}
		rSum=__vsumi(rArr, nLength);

		for(int j=0;j<nLength;j++){
			mEmmissionArr[i*nLength+j]=1.0*rArr[j]/rSum;
		}
	}

	free(rArr);

	// debug
	// {
	// 	float sum=0;
	// 	float *arr1=NULL;

	// 	arr1=__vnew(hmmObj->sLength, NULL);

	// 	sum=__vsum(piArr, hmmObj->sLength);
	// 	if(sum==1.0){
	// 		printf("init pi is sucess !!!\n");
	// 	}
		

	// 	__msum(mTransitionArr, hmmObj->sLength, hmmObj->sLength, 1, arr1);
	// 	sum=__vsum(arr1, hmmObj->sLength);
	// 	if(sum==hmmObj->sLength*1.0){
	// 		printf("init A matrix is sucess !!!\n");
	// 	}

	// 	__msum(mEmmissionArr, hmmObj->sLength, hmmObj->nLength, 1, arr1);
	// 	sum=__vsum(arr1, hmmObj->sLength);
	// 	if(sum==hmmObj->sLength*1.0){
	// 		printf("init B matrix is sucess !!!\n");
	// 	}
	// }
}

static float __hmmObj_predict(HMMObj hmmObj,int *oArr,int tLength,int type){
	float *piArr=NULL; // sLength
	float *mTransitionArr=NULL; // sLength*nLength
	float *mEmmissionArr=NULL; // sLength*nLength

	int sLength=0;
	int nLength=0;

	float *mArr=NULL;
	float pValue=0;

	piArr=hmmObj->piArr;
	mTransitionArr=hmmObj->mTransitionArr;
	mEmmissionArr=hmmObj->mEmmissionArr;

	sLength=hmmObj->sLength;
	nLength=hmmObj->nLength;

	mArr=__vnew(tLength*sLength, NULL);

	if(type==0){ // forward
		__forward(piArr,mTransitionArr,mEmmissionArr,
				sLength,nLength,
				mArr,oArr,
				tLength);

		pValue=__vsum(mArr+(tLength-1)*sLength, sLength);
	}
	else{
		__backward(piArr,mTransitionArr,mEmmissionArr,
				sLength,nLength,
				mArr,oArr,
				tLength);

		for(int i=0;i<sLength;i++){
			pValue+=mArr[i]*piArr[i]*mEmmissionArr[i*nLength+oArr[0]];
		}
	}
	
	free(mArr);
	return pValue;
}

static float __hmmObj_calGamma(HMMObj hmmObj,int i,int t){
	float gamma=0;
	float num=0;
	float den=0;

	int sLength=0;

	float *mAlphaArr=NULL; // forward/backward
	float *mBetaArr=NULL;

	sLength=hmmObj->sLength;

	mAlphaArr=hmmObj->mAlphaArr;
	mBetaArr=hmmObj->mBetaArr;

	num=mAlphaArr[t*sLength+i]*mBetaArr[t*sLength+i];
	for(int i=0;i<sLength;i++){
		den+=mAlphaArr[t*sLength+i]*mBetaArr[t*sLength+i];
	}

	gamma=num/(den);
	return gamma;
}

static float __hmmObj_calKsi(HMMObj hmmObj,int *oArr,int i,int j,int t){
	float ksi=0;
	float num=0;
	float den=0;

	int sLength=0;
	int nLength=0;

	float *mAlphaArr=NULL; // forward/backward
	float *mBetaArr=NULL;

	float *mTransitionArr=NULL;
	float *mEmmissionArr=NULL;

	sLength=hmmObj->sLength;
	nLength=hmmObj->nLength;

	mAlphaArr=hmmObj->mAlphaArr;
	mBetaArr=hmmObj->mBetaArr;

	mTransitionArr=hmmObj->mTransitionArr;
	mEmmissionArr=hmmObj->mEmmissionArr;

	num=mAlphaArr[t*sLength+i]*mTransitionArr[i*sLength+j]*
		mEmmissionArr[j*nLength+oArr[t+1]]*mBetaArr[(t+1)*sLength+j];

	for(int i=0;i<sLength;i++){
		for(int j=0;j<sLength;j++){
			den+=mAlphaArr[t*sLength+i]*mTransitionArr[i*sLength+j]*
					mEmmissionArr[j*nLength+oArr[t+1]]*mBetaArr[(t+1)*sLength+j];
		}
	}

	ksi=num/(den);
	return ksi;
}

// alpha matrix
static void __forward(float *piArr,float *mTransitionArr,float *mEmmissionArr,
					int sLength,int nLength,
					float *mAlphaArr,int *oArr,
					int tLength){
	float sum=0;

	// t0 column init
	for(int i=0;i<sLength;i++){
		mAlphaArr[i]=piArr[i]*mEmmissionArr[i*nLength+oArr[0]];
	}

	for(int i=1;i<tLength;i++){
		for(int j=0;j<sLength;j++){ // a[ij]*b[j(t+1)]*alpha[ti]
			sum=0;
			for(int k=0;k<sLength;k++){
				sum+=mAlphaArr[(i-1)*sLength+k]*
						mTransitionArr[k*sLength+j];
			}

			mAlphaArr[i*sLength+j]=sum*mEmmissionArr[j*nLength+oArr[i]];
		}
	}
}

// beta matrix
static void __backward(float *piArr,float *mTransitionArr,float *mEmmissionArr,
					int sLength,int nLength,
					float *mBetaArr,int *oArr,
					int tLength){
	float sum=0;

	// t-1 column init
	for(int i=0;i<sLength;i++){
		mBetaArr[(tLength-1)*sLength+i]=1;
	}

	for(int i=tLength-2;i>=0;i--){
		for(int j=0;j<sLength;j++){ // a[ij]*b[j(t+1)]*beta[(t+1)i]
			sum=0;
			for(int k=0;k<sLength;k++){
				sum+=mTransitionArr[j*sLength+k]*
						mEmmissionArr[k*nLength+oArr[i+1]]*
						mBetaArr[(i+1)*sLength+k];
			}

			mBetaArr[i*sLength+j]=sum;
		}
	}
}

// sum(arr)==1.0
static int __distribute(float *arr,int length){
	int index=0;
	float r1=0;

	time_t t;	
	srand((unsigned) time(&t));

	r1=rand()%10001/10000.0;
	for(int i=0;i<length;i++){
		if(r1<arr[i]+1e-4){
			index=i;
			break;
		}
		r1-=arr[i];
	}

	return index;
}	







