// 

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "viterbi.h"

static void __viterbi1(float *piArr,float *mAArr,float *mBArr,
						int sLength,int nLength,
						int *oArr,int tLength,
						float *mProbArr,int *mIndexArr);

static void __viterbi2(float *piArr,float *mAArr,float *mBArr,
						int sLength,int nLength,
						int *oArr,int tLength,
						float *mProbArr,int *mIndexArr);
/***
	λ=(π,A,B) sLength*1,sLength*sLength,sLength*nLength
	oArr B nLength index
	mProbArr tLength*sLength
	sArr A sLength index
****/
float viterbi(float *piArr,float *mAArr,float *mBArr,
			int sLength,int nLength,
			int *oArr,int tLength,int *isLog,
			int *sArr,float *mProbArr,int *mIndexArr){
	float prob=0;

	float *piArr1=NULL;
	float *mAArr1=NULL;
	float *mBArr1=NULL;

	int *oArr1=NULL;
	float *mProbArr1=NULL;

	int _isLog=0;

	if(isLog){
		_isLog=*isLog;
	}

	if(oArr){
		oArr1=oArr;
	}
	else{
		tLength=__varangei(0, nLength, 1, &oArr1);
	}

	if(mProbArr){
		mProbArr1=mProbArr;
	}
	else{
		mProbArr1=__vnew(tLength*sLength, NULL);
	}

	if(!_isLog){
		piArr1=piArr;
		mAArr1=mAArr;
		mBArr1=mBArr;
	}
	else{
		piArr1=__vnew(sLength, NULL);
		mAArr1=__vnew(sLength*sLength, NULL);
		mBArr1=__vnew(sLength*nLength, NULL);

		__vadd_value(piArr, 1e-16, sLength, piArr1);
		__vadd_value(mAArr, 1e-16, sLength*sLength, mAArr1);
		__vadd_value(mBArr, 1e-16, sLength*nLength, mBArr1);

		__vlog(piArr1, sLength, NULL);
		__vlog(mAArr1, sLength*sLength, NULL);
		__vlog(mBArr1, sLength*nLength, NULL);
	}

	// viterbi
	if(!_isLog){
		__viterbi1(piArr1,mAArr1,mBArr1,
				sLength,nLength,
				oArr1,tLength,
				mProbArr1,mIndexArr);
	}
	else{
		__viterbi2(piArr1,mAArr1,mBArr1,
				sLength,nLength,
				oArr1,tLength,
				mProbArr1,mIndexArr);
	}

	// find hidden states
	if(sArr){
		for(int i=0;i<tLength;i++){
			sArr[i]=__vmax(mProbArr1+i*sLength, sLength, NULL);
		}
	}

	prob=mProbArr1[(tLength-1)*sLength+sArr[tLength-1]];

	if(!oArr){
		free(oArr1);
	}

	if(!mProbArr){
		free(mProbArr1);
	}

	if(_isLog){
		free(piArr1);
		free(mAArr1);
		free(mBArr1);
	}

	return prob;
}

static void __viterbi1(float *piArr,float *mAArr,float *mBArr,
						int sLength,int nLength,
						int *oArr,int tLength,
						float *mProbArr,int *mIndexArr){
	int oIndex=0;
	float eValue=0;

	float *arr1=NULL;

	arr1=__vnew(sLength, NULL);

	// first max likelihood
	for(int j=0;j<sLength;j++){
		mProbArr[j]=piArr[j]*mBArr[j*nLength+oArr[0]];
	}

	// second 
	for(int i=1;i<tLength;i++){
		oIndex=oArr[i];
		for(int j=0;j<sLength;j++){
			float _max=0;
			int _index=0;

			for(int k=0;k<sLength;k++){
				arr1[k]=mProbArr[(i-1)*sLength+k]*mAArr[k*sLength+j];
			}

			_index=__vmax(arr1, sLength, &_max);
			eValue=mBArr[j*nLength+oIndex];

			mProbArr[i*sLength+j]=_max*eValue;
			if(mIndexArr){
				mIndexArr[i*sLength+j]=_index;
			}
		}	
	}

	free(arr1);
}

static void __viterbi2(float *piArr,float *mAArr,float *mBArr,
						int sLength,int nLength,
						int *oArr,int tLength,
						float *mProbArr,int *mIndexArr){
	int oIndex=0;
	float eValue=0;

	float *arr1=NULL;

	arr1=__vnew(sLength, NULL);

	// first max likelihood
	for(int j=0;j<sLength;j++){
		mProbArr[j]=piArr[j]+mBArr[j*nLength+oArr[0]];
	}

	// second 
	for(int i=1;i<tLength;i++){
		oIndex=oArr[i];
		for(int j=0;j<sLength;j++){
			float _max=0;
			int _index=0;

			for(int k=0;k<sLength;k++){
				arr1[k]=mProbArr[(i-1)*sLength+k]+mAArr[k*sLength+j];
			}

			_index=__vmax(arr1, sLength, &_max);
			eValue=mBArr[j*nLength+oIndex];

			mProbArr[i*sLength+j]=_max+eValue;
			if(mIndexArr){
				mIndexArr[i*sLength+j]=_index;
			}
		}	
	}

	free(arr1);
}







