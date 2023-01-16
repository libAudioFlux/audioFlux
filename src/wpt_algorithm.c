// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/conv_algorithm.h"

#include "filterbank/dwt_filterCoef.h"

#include "wpt_algorithm.h"

struct OpaqueWPT{
	ConvObj convObj;

	int num; // level <=radix2Exp-1
	int radix2Exp;
	int fftLength; // 1<<radix2Exp
	
	float *loDArr;
	float *hiDArr;

	int decLength;

	int *indexArr; // 1<<num
	float *mNodeArr; // (level+1)*dataLength

	// cache
	float *curDataArr;
	float *cAArr; //  approximation&detail
	float *cDArr;

	int cacheLength; // fftLength+2*decLength

};

static void __periodPadding(float *arr1,int length1,int filterLength,float *arr2);
static float *__getNodeOffset(float *mNodeArr,int nodeIndex,int dataLength,int *outLength);

static void __sortIndex(int *vArr1,int length,int type,int *iArr2);
static void __reassign(int *indexArr1,int length,int *indexArr2);

static void __calIndexArr(int num,int *indexArr);

/***
	num radix2Exp-1 <=radix2Exp-1
	waveletType 'sym4' 'db4'/'coif4'/'fk4'/'bior4.4'
****/
int wptObj_new(WPTObj *wptObj,int num,int radix2Exp,
			 WaveletDiscreteType *waveletType,int *t1,int *t2){
	int status=0;
	WPTObj wpt=NULL;

	ConvObj convObj=NULL;

	int fftLength=0;

	float *loDArr=NULL;
	float *hiDArr=NULL;

	int *indexArr=NULL; // 1<<num

	int decLength=0;
	int cacheLength=0;

	WaveletDiscreteType _waveletType=WaveletDiscrete_Sym;
	int _t1=4;
	int _t2=4;

	if(radix2Exp){
		if(radix2Exp<1||radix2Exp>30){
			status=-100;
			printf("radix2Exp is error!\n");
			return status;
		}
	}

	fftLength=1<<radix2Exp;

	if(waveletType){
		_waveletType=*waveletType;
	}

	if(num<1||num>radix2Exp-1){ 
		printf("num is error!\n");
		return -1;
	}

	if(t1){
		_t1=*t1;
	}

	if(t2){
		_t2=*t2;
	}

	wpt=*wptObj=(WPTObj )calloc(1, sizeof(struct OpaqueWPT ));

	convObj_new(&convObj);

	decLength=dwt_filterCoef(_waveletType,_t1,_t2,0,
							&loDArr,&hiDArr);

	cacheLength=fftLength+2*decLength;

	indexArr=__vnewi(1<<num, NULL);
	__calIndexArr(num, indexArr);

	wpt->convObj=convObj;

	wpt->num=num;
	wpt->radix2Exp=radix2Exp;
	wpt->fftLength=fftLength;

	wpt->loDArr=loDArr;
	wpt->hiDArr=hiDArr;

	wpt->decLength=decLength;

	wpt->indexArr=indexArr;
	wpt->mNodeArr=__vnew((num+1)*fftLength, NULL);

	wpt->curDataArr=__vnew(cacheLength, NULL);
	wpt->cAArr=__vnew(cacheLength, NULL);
	wpt->cDArr=__vnew(cacheLength, NULL);

	wpt->cacheLength=cacheLength;

	return status;
}

/***
	1. padding
	2. conv
	3. order
	4. reassign --> mDataArr
	coefArr dataLength;mDataArr 2^num*dataLength
****/
void wptObj_wpt(WPTObj wptObj,float *dataArr,float *coefArr,float *mDataArr){
	ConvObj convObj=NULL;

	int num=0; // level <=radix2Exp-1
	int radix2Exp=0;
	int fftLength=0; 

	float *loDArr=NULL;
	float *hiDArr=NULL;

	int decLength=0;

	int *indexArr=NULL;
	float *mNodeArr=NULL;

	float *curDataArr=NULL;

	float *cAArr=NULL; //  approximation&detail
	float *cDArr=NULL;

	int cacheLength=0;

	int downLength=0;
	int curDataLength=0;

	int count=0;
	int nodeIndex=1;

	indexArr=wptObj->indexArr;
	mNodeArr=wptObj->mNodeArr;

	convObj=wptObj->convObj;

	num=wptObj->num;
	radix2Exp=wptObj->radix2Exp;
	fftLength=wptObj->fftLength;

	loDArr=wptObj->loDArr;
	hiDArr=wptObj->hiDArr;

	decLength=wptObj->decLength;
	cacheLength=wptObj->cacheLength;

	curDataArr=wptObj->curDataArr;

	cAArr=wptObj->cAArr;
	cDArr=wptObj->cDArr;

	ConvModeType modeType=ConvMode_Valid;
	ConvMethodType methodType=ConvMethod_Auto;

	count=(1<<num)-1;

	memcpy(mNodeArr, dataArr, sizeof(float )*fftLength);

	// debug
	// {
	// 	printf("LoD is:\n");
	// 	__vdebug(loDArr, decLength, 1);
	// 	printf("\n");

	// 	printf("HiD is:\n");
	// 	__vdebug(hiDArr, decLength, 1);
	// 	printf("\n");
	// }

	for(int i=0;i<count;i++){
		float *p=0;

		p=__getNodeOffset(mNodeArr, i, fftLength, &downLength);

		// 1. padding
		memset(curDataArr, 0, sizeof(float )*cacheLength);
		__periodPadding(p, downLength,decLength, curDataArr);

		curDataLength=downLength+decLength;

		// 2. conv
		memset(cAArr, 0, sizeof(float )*cacheLength);
		convObj_conv(convObj,curDataArr,curDataLength,loDArr,decLength,
					&modeType,&methodType,
					cAArr);
		for(int i=0;i<downLength/2;i++){
			cAArr[i]=cAArr[i*2+1];
		}

		memset(cDArr, 0, sizeof(float )*cacheLength);
		convObj_conv(convObj,curDataArr,curDataLength,hiDArr,decLength,
					&modeType,&methodType,
					cDArr);
		for(int i=0;i<downLength/2;i++){
			cDArr[i]=cDArr[i*2+1];
		}

		// 3. order
		p=__getNodeOffset(mNodeArr, nodeIndex, fftLength, &downLength);
		memcpy(p, (i&&i%2==0?cDArr:cAArr), sizeof(float )*downLength);

		p=__getNodeOffset(mNodeArr, nodeIndex+1, fftLength, &downLength);
		memcpy(p, (i&&i%2==0?cAArr:cDArr), sizeof(float )*downLength);

		nodeIndex+=2;
	}

	for(int i=0;i<count+1;i++){
		// i -> indexArr[i]*downLength ???
		memcpy(coefArr+i*downLength, mNodeArr+(num*fftLength+i*downLength), sizeof(float )*downLength);
	}

	// 4. reassign
	if(mDataArr){
		int start=0;
		int end=0;

		int bLen=0;
		int kLen=0;

		bLen=downLength;
		for(int i=0;i<count+1;i++){
			start=i*downLength;
			end=(i+1)*downLength-1;

			kLen=fftLength/bLen;
			for(int k=0;k<kLen;k++){
				for(int j=k,l=start;j<fftLength;j+=kLen,l++){
					mDataArr[i*fftLength+j]=coefArr[l]; 
				}
			}
		}
	}
}

static void __calIndexArr(int num,int *indexArr){
	int base=0;
	int *arr=NULL;

	int length=0;

	length=1<<num;
	base=length-1;

	__varangei(base, base+length-1, 1, &arr);

	__reassign(arr, length, indexArr);

	free(arr);
}

static void __reassign(int *indexArr1,int length,int *indexArr2){
	int num=0;
	int *arr=NULL;

	int curIndex=1;
	int len=0;

	arr=__vnewi(length, NULL);
	num=floorf(log2f(length));
	for(int i=0;i<num;i++){
		len=1<<i;
		for(int j=0;j<len;j++){
			arr[curIndex]=(1<<(i+1))-1-arr[j];
			curIndex++;
		}
	}

	__sortIndex(arr,length,0,indexArr2);

	free(arr);
}

static void __sortIndex(int *vArr1,int length,int type,int *iArr2){
	int *arr=NULL;

	arr=__vnewi(length, NULL);
	memcpy(arr, vArr1, sizeof(int )*length);

	for(int i=0;i<length;i++){
		iArr2[i]=i;
	}

	for(int i=0;i<length;i++){
		for(int j=i+1;j<length;j++){
			int _value=0;
			int _index=0;

			if(!type){ // 升
				if(arr[i]>arr[j]){
					_value=arr[i];
					arr[i]=arr[j];
					arr[j]=_value;

					_index=iArr2[i];
					iArr2[i]=iArr2[j];
					iArr2[j]=_index;
				}
			}
			else{ // 降
				if(arr[i]<arr[j]){
					_value=arr[i];
					arr[i]=arr[j];
					arr[j]=_value;

					_index=iArr2[i];
					iArr2[i]=iArr2[j];
					iArr2[j]=_index;
				}
			}
		}
	}

	free(arr);
}

static float *__getNodeOffset(float *mNodeArr,int nodeIndex,int dataLength,int *outLength){
	float *p=NULL;
	int curLength=0;
	
	int base=0;

	if(nodeIndex){
		int j=0;

		base=floor(log2f(nodeIndex+1));
		j=nodeIndex+1-(1<<base);

		curLength=dataLength/(1<<base);
		p=mNodeArr+(base*dataLength+j*curLength);
	}
	else{
		curLength=dataLength;
		p=mNodeArr;
	}

	if(outLength){
		*outLength=curLength;
	}
	
	return p;
}

static void __periodPadding(float *arr1,int length1,int filterLength,float *arr2){
	int totalLen=0;
	int halfLen=0;

	int firstIndex=0;
	int lastIndex=0;
	int n=0;

	totalLen=length1+filterLength;
	halfLen=filterLength/2;
	if(length1>=halfLen){
		memcpy(arr2,arr1+(length1-halfLen),sizeof(float )*halfLen);
		memcpy(arr2+halfLen,arr1,sizeof(float )*length1);
		memcpy(arr2+(halfLen+length1),arr1,sizeof(float )*halfLen);
	}
	else{
		int curIndex=0;

		firstIndex=(length1-halfLen+1)%length1;
		if(firstIndex<0){
			firstIndex=length1+firstIndex;
		}
		
		if(firstIndex==0){
			firstIndex=length1-1;
		}
		else{
			firstIndex=firstIndex-1;
		}

		n=floorf((totalLen-(length1-firstIndex))*1.0/length1);

		lastIndex=totalLen-(length1-firstIndex)-n*length1-1;

		memcpy(arr2, arr1+firstIndex, sizeof(float )*(length1-firstIndex));
		curIndex=length1-firstIndex;
		for(int i=0;i<n;i++){
			memcpy(arr2+curIndex,arr1,sizeof(float )*length1);
			curIndex+=length1;
		}

		memcpy(arr2+curIndex, arr1, sizeof(float )*(lastIndex+1));
	}
}

void wptObj_free(WPTObj wptObj){
	ConvObj convObj=NULL;
	
	float *loDArr=NULL;
	float *hiDArr=NULL;

	int *indexArr=NULL;
	float *mNodeArr=NULL; 

	float *curDataArr=NULL;
	float *cAArr=NULL; //  approximation&detail
	float *cDArr=NULL;

	if(wptObj){
		convObj=wptObj->convObj;

		loDArr=wptObj->loDArr;
		hiDArr=wptObj->hiDArr;

		indexArr=wptObj->indexArr;
		mNodeArr=wptObj->mNodeArr;

		curDataArr=wptObj->curDataArr;
		cAArr=wptObj->cAArr;
		cDArr=wptObj->cDArr;

		convObj_free(convObj);

		free(loDArr);
		free(hiDArr);

		free(indexArr);
		free(mNodeArr);

		free(curDataArr);
		free(cAArr);
		free(cDArr);

		free(wptObj);
	}
}









