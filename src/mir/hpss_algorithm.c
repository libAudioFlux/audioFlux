// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_complex.h"
#include "../vector/flux_vectorOp.h"

#include "../util/flux_util.h"

#include "../stft_algorithm.h"

#include "hpss_algorithm.h"

struct OpaqueHPSS{
	STFTObj stftObj;

	int fftLength;
	int slideLength;

	int hOrder;
	int pOrder;

	int timeLength;

	float *mRealArr; // timeLength*fftLength
	float *mImageArr;

	float *mPhaseImageArr; // timeLength*(fftLength/2+1)
	float *mPhaseRealArr;

	float *mMagArr;

	float *mHArr; // filter&mask
	float *mPArr;

};

int hpssObj_new(HPSSObj *hpssObj,
				int radix2Exp,WindowType *windowType,int *slideLength,
				int *hOrder,int *pOrder){
	int status=0;

	int fftLength=0;
	int _slideLength=0;

	int _hOrder=21;
	int _pOrder=31;

	WindowType _windowType=Window_Hamm;

	HPSSObj hpss=NULL;
	STFTObj stftObj=NULL;

	hpss=*hpssObj=(HPSSObj )calloc(1,sizeof(struct OpaqueHPSS ));

	fftLength=1<<radix2Exp;
	if(slideLength){
		if(*slideLength>0){
			_slideLength=*slideLength;
		}
	}

	if(hOrder){
		if(*hOrder>0&&(*hOrder&1)){
			_hOrder=*hOrder;
		}
	}

	if(pOrder){
		if(*pOrder>0&&(*pOrder&1)){
			_pOrder=*pOrder;
		}
	}

	if(windowType){
		_windowType=*windowType;
	}

	_slideLength=fftLength/4;

	stftObj_new(&stftObj, radix2Exp, &_windowType, &_slideLength, NULL);
		
	hpss->stftObj=stftObj;

	hpss->fftLength=fftLength;
	hpss->slideLength=_slideLength;

	hpss->hOrder=_hOrder;
	hpss->pOrder=_pOrder;

	return status;
}

int hpssObj_calDataLength(HPSSObj hpssObj,int dataLength){
	int length=0;

	int fftLength=0;
	int slideLength=0;

	int timeLength=0;

	fftLength=hpssObj->fftLength;
	slideLength=hpssObj->slideLength;

	timeLength=stftObj_calTimeLength(hpssObj->stftObj, dataLength);
	length=(timeLength-1)*slideLength+fftLength;

	return length;
}

/***
	1. stft&phase
	2. filter&mask
	3. istft
****/
void hpssObj_hpss(HPSSObj hpssObj,float *dataArr,int dataLength,float *hArr,float *pArr){
	STFTObj stftObj=NULL;

	int fftLength=0;
	int slideLength=0;

	int hOrder=0;
	int pOrder=0;

	float *mRealArr=NULL;
	float *mImageArr=NULL;

	float *mPhaseImageArr=NULL;
	float *mPhaseRealArr=NULL;

	float *mMagArr=NULL;

	float *mHArr=NULL; // filter&mask
	float *mPArr=NULL;

	int timeLength=0;

	if(!hArr&&!pArr){
		return;
	}

	stftObj=hpssObj->stftObj;

	fftLength=hpssObj->fftLength;
	slideLength=hpssObj->slideLength;

	hOrder=hpssObj->hOrder;
	pOrder=hpssObj->pOrder;

	timeLength=stftObj_calTimeLength(stftObj, dataLength);
	if(timeLength>hpssObj->timeLength||
		hpssObj->timeLength>2*timeLength){

		free(hpssObj->mRealArr);
		free(hpssObj->mImageArr);

		free(hpssObj->mPhaseImageArr);
		free(hpssObj->mPhaseRealArr);

		free(hpssObj->mMagArr);

		free(hpssObj->mHArr);
		free(hpssObj->mPArr);

		hpssObj->mRealArr=__vnew(timeLength*fftLength, NULL);
		hpssObj->mImageArr=__vnew(timeLength*fftLength, NULL);

		hpssObj->mPhaseImageArr=__vnew(timeLength*(fftLength/2+1), NULL);
		hpssObj->mPhaseRealArr=__vnew(timeLength*(fftLength/2+1), NULL);

		hpssObj->mMagArr=__vnew(timeLength*(fftLength/2+1), NULL);

		hpssObj->mHArr=__vnew(timeLength*(fftLength/2+1), NULL);
		hpssObj->mPArr=__vnew(timeLength*(fftLength/2+1), NULL);

		hpssObj->timeLength=timeLength;
	}

	mRealArr=hpssObj->mRealArr;
	mImageArr=hpssObj->mImageArr;

	mPhaseImageArr=hpssObj->mPhaseImageArr;
	mPhaseRealArr=hpssObj->mPhaseRealArr;

	mMagArr=hpssObj->mMagArr;

	mHArr=hpssObj->mHArr;
	mPArr=hpssObj->mPArr;

	// 1. stft&phase
	stftObj_stft(stftObj, dataArr, dataLength, mRealArr, mImageArr);

	__mcsquare2(mRealArr, mImageArr, timeLength, fftLength, fftLength/2+1, mMagArr); // S^2
	__vsqrt(mMagArr, timeLength*(fftLength/2+1), NULL);

	for(int i=0;i<timeLength;i++){
		for(int j=0;j<fftLength/2+1;j++){
			float v1=0;

			float r1=0;
			float i1=0;

			r1=mRealArr[i*fftLength+j];
			i1=mImageArr[i*fftLength+j];

			v1=mMagArr[i*(fftLength/2+1)+j];
			if(v1<1e-16){
				v1=1e-16;
			}

			mPhaseRealArr[i*(fftLength/2+1)+j]=r1/v1;
			mPhaseImageArr[i*(fftLength/2+1)+j]=i1/v1;
		}
	}

	// {
	// 	printf("mag is :\n");
	// 	__mdebug(mMagArr, timeLength, fftLength/2+1, 1);
	// 	printf("\n");

	// 	printf("stft is :\n");
	// 	__mcdebug(mRealArr,mImageArr, timeLength, fftLength, 1);
	// 	printf("\n");

	// 	printf("phase is :\n");
	// 	__mcdebug(mPhaseRealArr,mPhaseImageArr, timeLength, fftLength/2+1, 1);
	// 	printf("\n");
	// }

	// 2. filter&mask
	__mmedianfilter(mMagArr, timeLength, fftLength/2+1,	0, hOrder, mHArr);
	__mmedianfilter(mMagArr, timeLength, fftLength/2+1, 1, pOrder, mPArr);

	// {
	// 	printf("mHArr1 is :\n");
	// 	__mdebug(mHArr, timeLength, fftLength/2+1, 1);
	// 	printf("\n");

	// 	printf("mPArr1 is :\n");
	// 	__mdebug(mPArr, timeLength, fftLength/2+1, 1);
	// 	printf("\n");
	// }

	for(int i=0;i<timeLength*(fftLength/2+1);i++){
		float h1=0;
		float p1=0;

		float v1=0;
		float v2=0;

		h1=mHArr[i]*mHArr[i];
		p1=mPArr[i]*mPArr[i];
		
		v1=h1+p1;
		if(v1<1e-16){
			v1=1e-16;
		}

		v2=mMagArr[i];

		mHArr[i]=h1/v1*v2;
		mPArr[i]=p1/v1*v2;
	}

	// {
	// 	printf("mHArr2 is :\n");
	// 	__mdebug(mHArr, timeLength, fftLength/2+1, 1);
	// 	printf("\n");

	// 	printf("mPArr2 is :\n");
	// 	__mdebug(mPArr, timeLength, fftLength/2+1, 1);
	// 	printf("\n");
	// }

	// 3. istft
	if(hArr){
		for(int i=0;i<timeLength;i++){
			for(int j=0;j<fftLength/2+1;j++){
				float v1=0;

				float r1=0;
				float i1=0;

				v1=mHArr[i*(fftLength/2+1)+j];

				r1=mPhaseRealArr[i*(fftLength/2+1)+j]*v1;
				i1=mPhaseImageArr[i*(fftLength/2+1)+j]*v1;

				mRealArr[i*fftLength+j]=r1;
				mImageArr[i*fftLength+j]=i1;

				if(j>0&&j<fftLength/2){
					mRealArr[i*fftLength+(fftLength-j)]=r1;
					mImageArr[i*fftLength+(fftLength-j)]=-i1;
				}
			}
		}
		stftObj_istft(stftObj, mRealArr, mImageArr, timeLength, 0, hArr);
	}
	
	if(pArr){
		for(int i=0;i<timeLength;i++){
			for(int j=0;j<fftLength/2+1;j++){
				float v1=0;

				float r1=0;
				float i1=0;

				v1=mPArr[i*(fftLength/2+1)+j];

				r1=mPhaseRealArr[i*(fftLength/2+1)+j]*v1;
				i1=mPhaseImageArr[i*(fftLength/2+1)+j]*v1;

				mRealArr[i*fftLength+j]=r1;
				mImageArr[i*fftLength+j]=i1;

				if(j>0&&j<fftLength/2){
					mRealArr[i*fftLength+(fftLength-j)]=r1;
					mImageArr[i*fftLength+(fftLength-j)]=-i1;
				}
			}
		}
		stftObj_istft(stftObj, mRealArr, mImageArr, timeLength, 0, pArr);
	}

	// {
	// 	if(hArr){
	// 		printf("hArr is :\n");
	// 		for(int i=0;i<dataLength;i++){
	// 			printf("%d: %f\n",i,hArr[i]);
	// 		}
	// 		printf("\n");
	// 	}

	// 	if(pArr){
	// 		printf("pArr is :\n");
	// 		for(int i=0;i<dataLength;i++){
	// 			printf("%d: %f\n",i,pArr[i]);
	// 		}
	// 		printf("\n");
	// 	}
	// }
}

void hpssObj_free(HPSSObj hpssObj){

	if(!hpssObj){
		return;
	}

	stftObj_free(hpssObj->stftObj);

	free(hpssObj->mRealArr);
	free(hpssObj->mImageArr);

	free(hpssObj->mPhaseImageArr);
	free(hpssObj->mPhaseRealArr);

	free(hpssObj->mMagArr);

	free(hpssObj->mHArr);
	free(hpssObj->mPArr);

	free(hpssObj);
}

void hpssObj_debug(HPSSObj hpssObj){

}









