// 

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "phase_vocoder.h"

/***
	mDataArr1 stft time*fftLength
	rate 0.5~2
	slideLength fftLength/4
****/
void phase_vocoder(float *mRealArr1,float *mImageArr1,int nLength,int mLength,int slideLength,float rate,
				float *mRealArr2,float *mImageArr2){
	float *timeArr=NULL;
	float *phiArr=NULL;
	float *phaseArr=NULL;

	float *arr1=NULL;
	float *arr2=NULL;

	float *arr3=NULL;
	float *arr4=NULL;

	int tLen=0;
	int fLen=0;

	fLen=mLength/2+1;
	tLen=ceilf(nLength/rate);
	__varange(0,nLength,rate,&timeArr);

	phiArr=__vlinspace(0, M_PI*slideLength, fLen, 0);

	phaseArr=__vnew(fLen, NULL);
	__vcangle(mRealArr1,mImageArr1,fLen,phaseArr);

	arr1=__vnew(fLen, NULL);
	arr2=__vnew(fLen, NULL);

	arr3=__vnew(fLen, NULL);
	arr4=__vnew(fLen, NULL);
	for(int i=0;i<tLen;i++){
		int k=0;
		float alpha=0;

		float *rArr1=NULL;
		float *iArr1=NULL;

		float *rArr2=NULL;
		float *iArr2=NULL;

		k=floorf(timeArr[i]);
		alpha=timeArr[i]-floorf(timeArr[i]);

		if(k<nLength){
			rArr1=mRealArr1+k*mLength;
			iArr1=mImageArr1+k*mLength;
		}
		else{
			memset(arr1,0,sizeof(float )*fLen);
			memset(arr2,0,sizeof(float )*fLen);
			rArr1=arr1;
			iArr1=arr2;
		}
		
		if(k+1<nLength){
			rArr2=mRealArr1+(k+1)*mLength;
			iArr2=mImageArr1+(k+1)*mLength;
		}
		else{
			memset(arr3,0,sizeof(float )*fLen);
			memset(arr4,0,sizeof(float )*fLen);
			rArr2=arr3;
			iArr2=arr4;
		}

		__vcabs(rArr1, iArr1, fLen, arr1);
		__vcabs(rArr2, iArr2, fLen, arr2);
		__vmul_value(arr1, (1-alpha), fLen, NULL);
		__vmul_value(arr2, alpha, fLen, NULL);

		// cal data
		for(int j=0;j<fLen;j++){
			mRealArr2[i*mLength+j]=(arr1[j]+arr2[j])*cosf(phaseArr[j]);
			mImageArr2[i*mLength+j]=(arr1[j]+arr2[j])*sinf(phaseArr[j]);
		}

		for(int j=fLen,l=mLength/2-1;j<mLength;j++,l--){
			mRealArr2[i*mLength+j]=mRealArr2[i*mLength+l];
			mImageArr2[i*mLength+j]=-mImageArr2[i*mLength+l];
		}

		// update phase
		__vcangle(rArr2,iArr2,fLen,arr1);
		__vcangle(rArr1,iArr1,fLen,arr2);
		for(int j=0;j<fLen;j++){
			float value=0;

			value=arr1[j]-arr2[j]-phiArr[j];
			value=value-2*M_PI*roundf(value/(2*M_PI));

			phaseArr[j]+=(phiArr[j]+value);
		}
	}

	free(timeArr);
	free(phiArr);
	free(phaseArr);

	free(arr1);
	free(arr2);
	free(arr3);
	free(arr4);
}









