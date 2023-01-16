// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/conv_algorithm.h"

#include "filterbank/dwt_filterCoef.h"

#include "swt_algorithm.h"

struct OpaqueSWT{
	ConvObj convObj;

	int num; 
	int fftLength; // fftLength%2^num=0
	
	float *loDArr;
	float *hiDArr;

	int decLength;

	// cache
	float *curDataArr;

	float *cAArr; //  approximation&detail
	float *cDArr;

	float *loDArr2; // 2^num*decLength
	float *hiDArr2;

	int cacheLength; // fftLength+2^num*decLength

};

// periodic&symmetric
static void __periodPadding(float *arr1,int length1,int filterLength,float *arr2);

/***
	num >=1
	fftLength fftLength%2^num=0
	waveletType 'sym4' 'db4'/'coif4'/'fk4'/'bior4.4'
****/
int swtObj_new(SWTObj *swtObj,int num,int length,
			 WaveletDiscreteType *waveletType,int *t1,int *t2){
	int status=0;
	SWTObj swt=NULL;

	ConvObj convObj=NULL;

	int fftLength=0;

	float *loDArr=NULL;
	float *hiDArr=NULL;

	int decLength=0;
	int cacheLength=0;

	WaveletDiscreteType _waveletType=WaveletDiscrete_Sym;
	int _t1=4;
	int _t2=4;

	fftLength=length;

	if(waveletType){
		_waveletType=*waveletType;
	}

	if(fftLength<(1<<num)||fftLength%(1<<num)){ 
		printf("length is error!\n");
		return -1;
	}

	if(t1){
		_t1=*t1;
	}

	if(t2){
		_t2=*t2;
	}

	swt=*swtObj=(SWTObj )calloc(1, sizeof(struct OpaqueSWT ));

	convObj_new(&convObj);

	decLength=dwt_filterCoef(_waveletType,_t1,_t2,0,
							&loDArr,&hiDArr);

	cacheLength=fftLength+(1<<num)*decLength+1;

	swt->convObj=convObj;

	swt->num=num;
	swt->fftLength=fftLength;

	swt->loDArr=loDArr;
	swt->hiDArr=hiDArr;

	swt->decLength=decLength;

	swt->curDataArr=__vnew(cacheLength, NULL);

	swt->cAArr=__vnew(cacheLength, NULL);
	swt->cDArr=__vnew(cacheLength, NULL);

	swt->loDArr2=__vnew((1<<num)*decLength, NULL);
	swt->hiDArr2=__vnew((1<<num)*decLength, NULL);

	swt->cacheLength=cacheLength;

	return status;
}

/***
	1. padding
	2. conv&&keep
	3. upsample
	num*dataLength,mDataArr1 app,mDataArr2 det
****/
void swtObj_swt(SWTObj swtObj,float *dataArr,float *mDataArr1,float *mDataArr2){
	ConvObj convObj=NULL;

	int num=0;
	int fftLength=0; 

	float *loDArr=NULL;
	float *hiDArr=NULL;

	int decLength=0;

	float *curDataArr=NULL;

	float *cAArr=NULL; //  approximation&detail
	float *cDArr=NULL;

	float *loDArr2=NULL; // 2^num*decLength
	float *hiDArr2=NULL;

	int cacheLength=0;

	int upLength=0;
	int curDataLength=0;

	convObj=swtObj->convObj;

	num=swtObj->num;
	fftLength=swtObj->fftLength;

	loDArr=swtObj->loDArr;
	hiDArr=swtObj->hiDArr;

	decLength=swtObj->decLength;
	cacheLength=swtObj->cacheLength;

	curDataArr=swtObj->curDataArr;

	cAArr=swtObj->cAArr;
	cDArr=swtObj->cDArr;

	loDArr2=swtObj->loDArr2;
	hiDArr2=swtObj->hiDArr2;

	ConvModeType modeType=ConvMode_Full;
	ConvMethodType methodType=ConvMethod_Auto;

	upLength=decLength;
	memcpy(mDataArr1, dataArr, sizeof(float )*fftLength);

	memcpy(loDArr2, loDArr, sizeof(float )*upLength);
	memcpy(hiDArr2, hiDArr, sizeof(float )*upLength);

	// debug
	// {
	// 	printf("LoD is:\n");
	// 	__vdebug(loDArr, decLength, 1);
	// 	printf("\n");

	// 	printf("HiD is:\n");
	// 	__vdebug(hiDArr, decLength, 1);
	// 	printf("\n");
	// }

	for(int i=0;i<num;i++){

		// 1. padding
		memset(curDataArr, 0, sizeof(float )*cacheLength);
		__periodPadding(mDataArr1+(i?i-1:0)*fftLength, fftLength,upLength, curDataArr);

		curDataLength=fftLength+upLength;

		// 2. conv&keep
		memset(cAArr, 0, sizeof(float )*cacheLength);
		convObj_conv(convObj,curDataArr,curDataLength,loDArr2,upLength,
					&modeType,&methodType,
					cAArr);
		memcpy(mDataArr1+i*fftLength, cAArr+upLength, sizeof(float )*fftLength);

		memset(cDArr, 0, sizeof(float )*cacheLength);
		convObj_conv(convObj,curDataArr,curDataLength,hiDArr2,upLength,
					&modeType,&methodType,
					cDArr);
		memcpy(mDataArr2+i*fftLength, cDArr+upLength, sizeof(float )*fftLength);
		
		// 3. upsample
		memset(loDArr2, 0, sizeof(float )*upLength*2);
		memset(hiDArr2, 0, sizeof(float )*upLength*2);
		for(int j=0;j<decLength;j++){
			loDArr2[j*(1<<(i+1))]=loDArr[j];
			hiDArr2[j*(1<<(i+1))]=hiDArr[j];
		}

		upLength*=2;
	}
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

void swtObj_free(SWTObj swtObj){
	ConvObj convObj=NULL;
	
	float *loDArr=NULL;
	float *hiDArr=NULL;

	float *curDataArr=NULL;

	float *cAArr=NULL; //  approximation&detail
	float *cDArr=NULL;

	float *loDArr2=NULL; // 2^num*decLength
	float *hiDArr2=NULL;

	if(swtObj){
		convObj=swtObj->convObj;

		loDArr=swtObj->loDArr;
		hiDArr=swtObj->hiDArr;

		curDataArr=swtObj->curDataArr;

		cAArr=swtObj->cAArr;
		cDArr=swtObj->cDArr;

		loDArr2=swtObj->loDArr2;
		hiDArr2=swtObj->hiDArr2;

		convObj_free(convObj);

		free(loDArr);
		free(hiDArr);

		free(curDataArr);

		free(cAArr);
		free(cDArr);

		free(loDArr2);
		free(hiDArr2);

		free(swtObj);
	}
}









