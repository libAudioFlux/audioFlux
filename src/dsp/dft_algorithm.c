// clang -g -c 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_complex.h"

#include "dft_algorithm.h"

struct OpaqueDFT{
	int length; // data length

	float *realArr; // dft input r,i cache
	float *imageArr;

	double *mCosArr; // length*length cache
	double *mSinArr;

};

static void __dot(double *mRealArr1,double *mImageArr1,
				float *mRealArr2,float *mImageArr2,
				int nLength1,int mLength1,
				int nLength2,int mLength2,
				float *mRealArr3,float *mImageArr3);

int dftObj_new(DFTObj *dftObj,int length){
	int status=0;
	DFTObj dft=NULL;

	float *realArr=NULL; // dft input r,i cache
	float *imageArr=NULL;

	double *mCosArr=NULL; // fft w cache
	double *mSinArr=NULL;

	double value=0;

	dft=*dftObj=(DFTObj )calloc(1, sizeof(struct OpaqueDFT ));

	realArr=__vnew(length, NULL);
	imageArr=__vnew(length, NULL);

	mCosArr=(double *)calloc(length*length, sizeof(double ));
	mSinArr=(double *)calloc(length*length, sizeof(double ));
	for(int i=0;i<length;i++){
		for(int j=0;j<length;j++){
			value=2*M_PI*i*j/length;
			mCosArr[i*length+j]=cos(value);
			mSinArr[i*length+j]=-sin(value);
		}
	}

	dft->length=length;

	dft->realArr=realArr;
	dft->imageArr=imageArr;

	dft->mCosArr=mCosArr;
	dft->mSinArr=mSinArr;

	return status;
}

// Y=W*X nLen*nLen@nLen*1 => nLen*1
void dftObj_dft(DFTObj dftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2){
	int length=0;

	float *realArr=NULL;
	float *imageArr=NULL;

	double *mCosArr=NULL; // length*length cache
	double *mSinArr=NULL;

	length=dftObj->length;

	realArr=dftObj->realArr;
	imageArr=dftObj->imageArr;

	// ??? 
	if(realArr1){
		memcpy(realArr, realArr1, sizeof(float )*length);
	}
	else{
		memset(realArr, 0, sizeof(float )*length);
	}

	if(imageArr1){
		memcpy(imageArr, imageArr1, sizeof(float )*length);
	}
	else{
		memset(imageArr, 0, sizeof(float )*length);
	}

	mCosArr=dftObj->mCosArr;
	mSinArr=dftObj->mSinArr;

	__dot(mCosArr,mSinArr,
			realArr,imageArr,
			length,length,
			length,1,
			realArr2,imageArr2);
}

void dftObj_idft(DFTObj dftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2){
	int length=0;

	float *realArr=NULL;
	float *imageArr=NULL;

	double *mCosArr=NULL; // length*length cache
	double *mSinArr=NULL;

	length=dftObj->length;

	realArr=dftObj->realArr;
	imageArr=dftObj->imageArr;

	// ??? 
	if(realArr1){
		memcpy(realArr, realArr1, sizeof(float )*length);
	}
	else{
		memset(realArr, 0, sizeof(float )*length);
	}

	if(imageArr1){
		memcpy(imageArr, imageArr1, sizeof(float )*length);
	}
	else{
		memset(imageArr, 0, sizeof(float )*length);
	}

	for(int i=0;i<length;i++){
		imageArr[i]=-imageArr[i];
	}

	mCosArr=dftObj->mCosArr;
	mSinArr=dftObj->mSinArr;

	__dot(mCosArr,mSinArr,
			realArr,imageArr,
			length,length,
			length,1,
			realArr2,imageArr2);

	for(int i=0;i<length;i++){
		realArr2[i]/=length;
		imageArr2[i]/=-length;
	}
}

void dftObj_free(DFTObj dftObj){
	float *realArr=NULL; // dft input r,i cache
	float *imageArr=NULL;

	double *mCosArr=NULL; // length*length cache
	double *mSinArr=NULL;

	if(!dftObj){
		return;
	}

	realArr=dftObj->realArr;
	imageArr=dftObj->imageArr;

	mCosArr=dftObj->mCosArr;
	mSinArr=dftObj->mSinArr;

	free(realArr);
	free(imageArr);

	free(mCosArr);
	free(mSinArr);

	free(dftObj);
}

static void __dot(double *mRealArr1,double *mImageArr1,
				float *mRealArr2,float *mImageArr2,
				int nLength1,int mLength1,
				int nLength2,int mLength2,
				float *mRealArr3,float *mImageArr3){
	if(mLength1!=nLength2){
		return;
	}

	for(int i=0;i<nLength1;i++){
		for(int j=0;j<mLength2;j++){
			double _value1=0;
			double _value2=0;

			double r1=0;
			double r2=0;

			double i1=0;
			double i2=0;

			for(int k=0;k<mLength1;k++){ // arr1.row*arr2.col
				r1=mRealArr1[i*mLength1+k];
				r2=mRealArr2[j+k*mLength2];

				i1=mImageArr1[i*mLength1+k];
				i2=mImageArr2[j+k*mLength2];

				_value1+=(r1*r2-i1*i2);
				_value2+=(i1*r2+r1*i2);
			}

			mRealArr3[i*mLength2+j]=_value1;
			mImageArr3[i*mLength2+j]=_value2;
		}
	}
}







