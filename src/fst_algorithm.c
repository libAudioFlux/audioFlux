// clang 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"


#include "dsp/flux_window.h"
#include "dsp/fft_algorithm.h"

#include "fst_algorithm.h"

struct OpaqueFST{
	FFTObj fftObj;
	FFTObj *fftObjArr; // radix2Exp-2

	int fftLength;
	int radix2Exp;

	float norm; // 1/sqrt(fftLength) 

	int *lenArr; // radix2Exp*2

	int *mIndexArr; // reassign fre*time -> fftLength*fftLength

	float *realArr1; // fftLength
	float *imageArr1;

	float *realArr2; 
	float *imageArr2;

	float *curDataArr; // ifftshift cache

};

static void _fstObj_initFFT(FSTObj fstObj,int radix2Exp);
static void _fstObj_initPartition(FSTObj fstObj,int length);
static void _fstObj_initReassign(FSTObj fstObj,int radix2Exp);

static void __mget_fast(float *mRealArr1,float *mImageArr1, int *mIndexArr,
				int length,int start,int end,
				float *mRealArr3,float *mImageArr3);

int fstObj_new(FSTObj *fstObj,int radix2Exp){
	int status=0;

	FSTObj fst=NULL;

	FFTObj fftObj=NULL;
	int fftLength=0;

	float norm=0;

	float *realArr1=NULL; // fftLength
	float *imageArr1=NULL;

	float *realArr2=NULL; 
	float *imageArr2=NULL;

	float *curDataArr=NULL;

	if(radix2Exp<3){ // dataLength>=8
		return -1;
	}

	fst=*fstObj=(FSTObj )calloc(1, sizeof(struct OpaqueFST ));
	
	fftLength=(1<<radix2Exp);
	fftObj_new(&fftObj, radix2Exp);

	realArr1=__vnew(fftLength, NULL);
	imageArr1=__vnew(fftLength, NULL);

	realArr2=__vnew(fftLength, NULL);
	imageArr2=__vnew(fftLength, NULL);

	curDataArr=__vnew(fftLength, NULL);
	
	_fstObj_initFFT(fst,radix2Exp);
	_fstObj_initPartition(fst,radix2Exp*2);
	_fstObj_initReassign(fst,radix2Exp);

	norm=1/sqrtf(fftLength);

	fst->fftObj=fftObj;

	fst->fftLength=fftLength;
	fst->radix2Exp=radix2Exp;

	fst->norm=norm;

	fst->realArr1=realArr1;
	fst->imageArr1=imageArr1;

	fst->realArr2=realArr2;
	fst->imageArr2=imageArr2;

	fst->curDataArr=curDataArr;

	return status;
}

/***
	1. fft -> out 
	2. 分段ifft
	3. 根据mIndex排列矩阵
****/
void fstObj_fst(FSTObj fstObj,float *dataArr,int minIndex,int maxIndex,float *mRealArr,float *mImageArr){
	FFTObj fftObj=NULL;
	FFTObj *fftObjArr=NULL; 

	int fftLength=0;

	float norm=0;

	int *lenArr=NULL; // radix2Exp*2
	int radix2Exp=0;

	int *mIndexArr=NULL; // reassign fre*time -> fftLength*fftLength

	float *realArr1=NULL; // fftLength
	float *imageArr1=NULL;

	float *realArr2=NULL; // fftLength
	float *imageArr2=NULL;

	float *curDataArr=NULL;

	int nLen=0;
	int mLen=0;

	int index=0;
	int len=0;

	fftObj=fstObj->fftObj;
	fftObjArr=fstObj->fftObjArr;

	fftLength=fstObj->fftLength;

	norm=fstObj->norm;

	lenArr=fstObj->lenArr;
	radix2Exp=fstObj->radix2Exp;

	mIndexArr=fstObj->mIndexArr;

	realArr1=fstObj->realArr1;
	imageArr1=fstObj->imageArr1;

	realArr2=fstObj->realArr2;
	imageArr2=fstObj->imageArr2;

	curDataArr=fstObj->curDataArr;

	if(minIndex<0){
		minIndex=0;
	}

	if(maxIndex>fftLength/2){
		maxIndex=fftLength/2;
	}

	if(minIndex>maxIndex){
		minIndex=0;
		maxIndex=fftLength/2;
	}

	// 1. fft ->out
	// ifftshift
	memcpy(curDataArr, dataArr+fftLength/2, sizeof(float )*fftLength/2);
	memcpy(curDataArr+fftLength/2, dataArr, sizeof(float )*fftLength/2);

	// fft
	fftObj_fft(fftObj, curDataArr, NULL, realArr1, imageArr1);

	// fftshift
	memcpy(realArr2, realArr1, sizeof(float )*fftLength);
	memcpy(imageArr2, imageArr1, sizeof(float )*fftLength);

	memcpy(realArr1, realArr2+fftLength/2, sizeof(float )*fftLength/2);
	memcpy(realArr1+fftLength/2, realArr2, sizeof(float )*fftLength/2);
	
	memcpy(imageArr1, imageArr2+fftLength/2, sizeof(float )*fftLength/2);
	memcpy(imageArr1+fftLength/2, imageArr2, sizeof(float )*fftLength/2);

	// norm
	__vmul_value(realArr1, norm, fftLength, NULL);
	__vmul_value(imageArr1, norm, fftLength, NULL);

	// 2. 分段ifft
	index=1;
	for(int i=1,j=0;i<radix2Exp-1;i++,j++){
		float _norm=0;

		len=lenArr[i];
		_norm=sqrtf(len);

		// ifftshift
		memcpy(realArr2, realArr1+index, sizeof(float )*len);
		memcpy(imageArr2, imageArr1+index, sizeof(float )*len);

		memcpy(realArr1+index, realArr2+len/2, sizeof(float )*len/2);
		memcpy(realArr1+(index+len/2), realArr2, sizeof(float )*len/2);

		memcpy(imageArr1+index, imageArr2+len/2, sizeof(float )*len/2);
		memcpy(imageArr1+(index+len/2), imageArr2, sizeof(float )*len/2);

		// ifft
		fftObj_ifft(fftObjArr[j], realArr1+index, imageArr1+index, realArr2, imageArr2);

		// norm
		__vmul_value(realArr2, _norm, len, NULL);
		__vmul_value(imageArr2, _norm, len, NULL);

		// fftshift
		memcpy(realArr1+index, realArr2+len/2, sizeof(float )*len/2);
		memcpy(realArr1+(index+len/2), realArr2, sizeof(float )*len/2);

		memcpy(imageArr1+index, imageArr2+len/2, sizeof(float )*len/2);
		memcpy(imageArr1+(index+len/2), imageArr2, sizeof(float )*len/2);

		index+=len;
	}

	index+=3;
	for(int i=radix2Exp+2,j=radix2Exp-3;i<2*radix2Exp;i++,j--){
		float _norm=0;

		len=lenArr[i];
		_norm=sqrtf(len);

		// ifftshift
		memcpy(realArr2, realArr1+index, sizeof(float )*len);
		memcpy(imageArr2, imageArr1+index, sizeof(float )*len);

		memcpy(realArr1+index, realArr2+len/2, sizeof(float )*len/2);
		memcpy(realArr1+(index+len/2), realArr2, sizeof(float )*len/2);

		memcpy(imageArr1+index, imageArr2+len/2, sizeof(float )*len/2);
		memcpy(imageArr1+(index+len/2), imageArr2, sizeof(float )*len/2);

		// ifft
		fftObj_ifft(fftObjArr[j], realArr1+index, imageArr1+index, realArr2, imageArr2);

		// norm
		__vmul_value(realArr2, _norm, len, NULL);
		__vmul_value(imageArr2, _norm, len, NULL);

		// fftshift
		memcpy(realArr1+index, realArr2+len/2, sizeof(float )*len/2);
		memcpy(realArr1+(index+len/2), realArr2, sizeof(float )*len/2);

		memcpy(imageArr1+index, imageArr2+len/2, sizeof(float )*len/2);
		memcpy(imageArr1+(index+len/2), imageArr2, sizeof(float )*len/2);

		index+=len;
	}

	// 3. 根据mIndex排列矩阵 fre*time -> (fftLength/2+1)*fftLength
	nLen=fftLength;
	mLen=fftLength;
	for(int i=fftLength/2-minIndex,k=0;i>=fftLength/2-maxIndex;i--,k++){
		for(int j=0;j<mLen;j++){
			int _index=0;

			_index=mIndexArr[(long long )i*mLen+j];
			mRealArr[(long long )k*mLen+j]=realArr1[_index];
			mImageArr[(long long )k*mLen+j]=imageArr1[_index];
		}
	}

	// __mget_fast(realArr1,imageArr1,mIndexArr,
	// 			fftLength,minIndex,maxIndex,
	// 			mRealArr,mImageArr);
}

static void _fstObj_initFFT(FSTObj fstObj,int radix2Exp){
	FFTObj *fftObjArr=NULL;

	fftObjArr=(FFTObj *)calloc(radix2Exp-2, sizeof(FFTObj ));
	for(int i=radix2Exp-2,j=0;i>=1;i--,j++){
		fftObj_new(fftObjArr+j, i);
	}

	fstObj->fftObjArr=fftObjArr;
}

static void _fstObj_initPartition(FSTObj fstObj,int length){
	int *lenArr=NULL;

	int *arr1=NULL;
	int value=0;

	value=length/2-1;
	__varangei(0, value, 1, &arr1);

	lenArr=__vnewi(length, NULL);

	lenArr[0]=1;
	lenArr[length/2-1]=1;
	lenArr[length/2]=1;
	for(int i=1,j=value-1;i<length/2-1;i++,j--){
		lenArr[i]=powf(2, arr1[j]);
	}

	for(int i=length/2+1,j=0;i<length;i++,j++){
		lenArr[i]=powf(2, arr1[j]);
	}

	free(arr1);
	fstObj->lenArr=lenArr;
}

static void _fstObj_initReassign(FSTObj fstObj,int radix2Exp){
	int fftLength=0;

	int *lenArr=NULL; 
	int *mIndexArr=NULL; 

	int value=0;

	fftLength=(1<<radix2Exp);

	lenArr=fstObj->lenArr;
	// mIndexArr=__vnewi(fftLength*fftLength, 0);
	/***
		1. optimize memory size
		2. long long resolve numeric overflow
	****/
	mIndexArr=(int *)calloc((long long )(fftLength/2+1)*fftLength, sizeof(int ));
	for(int i=0;i<radix2Exp*2;i++){
		int len1=0;
		int len2=0;

		int index1=0;
		int index2=0;

		len1=lenArr[i];
		len2=fftLength/len1;
		for(int j=0;j<len1;j++){
			index1=fftLength-__vsumi(lenArr,i+1);
			index2=len2*j;

			for(int k=index1;k<index1+len1;k++){
				if(k<fftLength/2+1){ // optimize >=fftLength/2+1 
					for(int l=index2;l<index2+len2;l++){ 
						mIndexArr[(long long)k*fftLength+l]=value;
					}
				}
			}

			value++;
		}
	}

	fstObj->mIndexArr=mIndexArr;
}

/***
	block per 200 for CPU cache
****/
static void __mget_fast(float *mRealArr1,float *mImageArr1, int *mIndexArr,
						int length,int start,int end,
						float *mRealArr3,float *mImageArr3){
	int nLen=0;

	int block=400;
	
	int num=0;
	int mod=0;

	int k=0;

	nLen=end-start+1;
	num=nLen/block;
	mod=nLen%block;

	for(int n=0;n<num;n++){
		for(int i=length/2-start-n*block;i>length/2-start-(n+1)*block;i--){
			for(int j=0;j<length;j++){
				int _index=0;

				_index=mIndexArr[i*length+j];
				mRealArr3[k*length+j]=mRealArr1[_index];
				mImageArr3[k*length+j]=mImageArr1[_index];
			}

			k++;
		}
	}

	if(mod){
		for(int i=length/2-end+mod-1;i>=length/2-end;i--){
			for(int j=0;j<length;j++){
				int _index=0;

				_index=mIndexArr[i*length+j];
				mRealArr3[k*length+j]=mRealArr1[_index];
				mImageArr3[k*length+j]=mImageArr1[_index];
			}

			k++;
		}
	}
}

void fstObj_free(FSTObj fstObj){
	FFTObj fftObj=NULL;
	FFTObj *fftObjArr=NULL; // radix2Exp-2

	int fftLength=0;

	int *lenArr=NULL; // radix2Exp*2
	int *mIndexArr=NULL; // reassign fre*time -> fftLength*fftLength

	float *realArr1=NULL; // fftLength
	float *imageArr1=NULL;

	float *realArr2=NULL; 
	float *imageArr2=NULL;

	float *curDataArr=NULL;

	if(!fstObj){
		return;
	}

	fftObj=fstObj->fftObj;
	fftObjArr=fstObj->fftObjArr;

	fftLength=fstObj->fftLength;

	lenArr=fstObj->lenArr;
	mIndexArr=fstObj->mIndexArr;

	realArr1=fstObj->realArr1;
	imageArr1=fstObj->imageArr1;

	realArr2=fstObj->realArr2;
	imageArr2=fstObj->imageArr2;

	curDataArr=fstObj->curDataArr;

	fftObj_free(fftObj);

	for(int i=0;i<fstObj->radix2Exp-2;i++){
		fftObj_free(fftObjArr[i]);
	}
	free(fftObjArr);

	free(lenArr);
	free(mIndexArr);

	free(realArr1);
	free(imageArr1);

	free(realArr2);
	free(imageArr2);

	free(curDataArr);

	free(fstObj);
}	











