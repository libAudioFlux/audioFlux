// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/conv_algorithm.h"

#include "filterbank/dwt_filterCoef.h"

#include "dwt_algorithm.h"

struct OpaqueDWT{
	ConvObj convObj;

	int num; // level <=radix2Exp-1
	int radix2Exp;
	int fftLength; // 1<<radix2Exp
	
	float *loDArr;
	float *hiDArr;

	int decLength;

	// cache
	float *curDataArr;
	float *cAArr; //  approximation&detail
	float *cDArr;

	int cacheLength; // fftLength+2*decLength

	// obsolete ???
	int samplate; 

	int *binBandArr;
	float *freBandArr;

};

// periodic&symmetric
static void __periodPadding(float *arr1,int length1,int filterLength,float *arr2);

static void __dwtObj_init(DWTObj dwtObj,WaveletDiscreteType waveletType,int t1,int t2);

/***
	num radix2Exp-1 <=radix2Exp-1
	samplate 32000
	waveletType 'sym4'
****/
int dwtObj_new(DWTObj *dwtObj,int num,int radix2Exp,
			 WaveletDiscreteType *waveletType,int *t1,int *t2){
	int status=0;
	DWTObj dwt=NULL;

	ConvObj convObj=NULL;

	int fftLength=0;
	int _samplate=32000;
	
	int *binBandArr=NULL;
	float *freBandArr=NULL;

	float *loDArr=NULL;
	float *hiDArr=NULL;

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

	dwt=*dwtObj=(DWTObj )calloc(1, sizeof(struct OpaqueDWT ));

	binBandArr=__vnewi(num, NULL);
	freBandArr=__vnew(num, NULL);
	for(int i=0;i<num;i++){
		binBandArr[i]=(1<<(i+1));
		freBandArr[i]=1.0*_samplate/fftLength*(1<<(i+1));
	}

	convObj_new(&convObj);

	decLength=dwt_filterCoef(_waveletType,_t1,_t2,0,
							&loDArr,&hiDArr);

	cacheLength=fftLength+2*decLength;

	dwt->convObj=convObj;

	dwt->num=num;
	dwt->radix2Exp=radix2Exp;
	dwt->fftLength=fftLength;

	dwt->samplate=_samplate;

	dwt->binBandArr=binBandArr;
	dwt->freBandArr=freBandArr;

	dwt->loDArr=loDArr;
	dwt->hiDArr=hiDArr;

	dwt->decLength=decLength;

	dwt->curDataArr=__vnew(cacheLength, NULL);
	dwt->cAArr=__vnew(cacheLength, NULL);
	dwt->cDArr=__vnew(cacheLength, NULL);

	dwt->cacheLength=cacheLength;

	return status;
}

static void __dwtObj_init(DWTObj dwtObj,WaveletDiscreteType waveletType,int t1,int t2){
	float *loDArr=NULL;
	float *hiDArr=NULL;

	int decLength=0;

	decLength=dwt_filterCoef(waveletType,t1,t2,0,
							&loDArr,&hiDArr);

	dwtObj->loDArr=loDArr;
	dwtObj->hiDArr=hiDArr;

	dwtObj->decLength=decLength;
}

float *dwtObj_getFreBandArr(DWTObj dwtObj){

	return dwtObj->freBandArr;
}

int *dwtObj_getBinBandArr(DWTObj dwtObj){

	return dwtObj->binBandArr;
}

/***
	1. padding
	2. conv
	3. split --> coefArr
	4. reassign --> mDataArr
	coefArr dataLength;mDataArr num*dataLength
****/
void dwtObj_dwt(DWTObj dwtObj,float *dataArr,float *coefArr,float *mDataArr){
	ConvObj convObj=NULL;

	int num=0; // level <=radix2Exp-1
	int radix2Exp=0;
	int fftLength=0; 

	float *loDArr=NULL;
	float *hiDArr=NULL;

	int decLength=0;

	float *curDataArr=NULL;

	float *cAArr=NULL; //  approximation & detail
	float *cDArr=NULL;

	int cacheLength=0;

	int downLength=0;
	int curDataLength=0;

	int cLength=0;

	convObj=dwtObj->convObj;

	num=dwtObj->num;
	radix2Exp=dwtObj->radix2Exp;
	fftLength=dwtObj->fftLength;

	loDArr=dwtObj->loDArr;
	hiDArr=dwtObj->hiDArr;

	decLength=dwtObj->decLength;
	cacheLength=dwtObj->cacheLength;

	curDataArr=dwtObj->curDataArr;

	cAArr=dwtObj->cAArr;
	cDArr=dwtObj->cDArr;

	ConvModeType modeType=ConvMode_Valid;
	ConvMethodType methodType=ConvMethod_Auto;

	downLength=fftLength;
	memcpy(cAArr, dataArr, sizeof(float )*downLength);

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
		__periodPadding(cAArr, downLength,decLength, curDataArr);

		curDataLength=downLength+decLength;
		downLength/=2;
		cLength+=downLength;

		// 2. conv
		memset(cAArr, 0, sizeof(float )*cacheLength);
		convObj_conv(convObj,curDataArr,curDataLength,loDArr,decLength,
					&modeType,&methodType,
					cAArr);
		for(int i=0;i<downLength;i++){
			cAArr[i]=cAArr[i*2+1];
		}

		memset(cDArr, 0, sizeof(float )*cacheLength);
		convObj_conv(convObj,curDataArr,curDataLength,hiDArr,decLength,
					&modeType,&methodType,
					cDArr);
		for(int i=0;i<downLength;i++){
			cDArr[i]=cDArr[i*2+1];
		}

		// 3. split
		memcpy(coefArr+(fftLength-cLength), cDArr, sizeof(float )*downLength);

		// debug
		// {
		// 	printf("datArr is :\n");
		// 	__vdebug(curDataArr, curDataLength, 1);
		// 	printf("\n");

		// 	printf("cAArr is :\n");
		// 	__vdebug(cAArr, downLength, 1);
		// 	printf("\n");

		// 	printf("cDArr is :\n");
		// 	__vdebug(cDArr, downLength, 1);
		// 	printf("\n");
		// }
	}
	memcpy(coefArr, cAArr, sizeof(float )*downLength);

	// 4. reassign
	if(mDataArr){
		int start=0;
		int end=0;

		int bLen=0;
		int kLen=0;

		for(int i=num;i>=1;i--){
			start=(1<<i);
			end=(1<<(i+1))-1;

			bLen=end-start+1;
			kLen=fftLength/bLen;

			for(int k=0;k<kLen;k++){
				for(int j=k,l=start;j<fftLength;j+=kLen,l++){
					mDataArr[(i-1)*fftLength+j]=coefArr[l]; // (num-i) ???
				}
			}
		}
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

void dwtObj_free(DWTObj dwtObj){
	ConvObj convObj=NULL;

	int *binBandArr=NULL;
	float *freBandArr=NULL;
	
	float *loDArr=NULL;
	float *hiDArr=NULL;

	float *curDataArr=NULL;
	float *cAArr=NULL; //  approximation&detail
	float *cDArr=NULL;

	if(dwtObj){
		convObj=dwtObj->convObj;

		binBandArr=dwtObj->binBandArr;
		freBandArr=dwtObj->freBandArr;

		loDArr=dwtObj->loDArr;
		hiDArr=dwtObj->hiDArr;

		curDataArr=dwtObj->curDataArr;
		cAArr=dwtObj->cAArr;
		cDArr=dwtObj->cDArr;

		convObj_free(convObj);

		free(binBandArr);
		free(freBandArr);

		free(loDArr);
		free(hiDArr);

		free(curDataArr);
		free(cAArr);
		free(cDArr);

		free(dwtObj);
	}
}









