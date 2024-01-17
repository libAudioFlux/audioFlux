// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "../util/flux_util.h"

#include "flux_window.h"
#include "filterDesign_fir.h"
#include "resample_algorithm.h"

struct OpaqueResample{
	int isContinue;
	int isScale;

	WindowType winType;
	float value; // kaiser/gauss
	float rollOff;

	int zeroNum;
	int bitLength; // 1<<nbit
	int interpLength; // zero*bitLength+1

	float *interpArr; // ratio<1 *ration
	float *interpDeltaArr;

	float ratio; // targetRate/sourceRate

	int p; // up
	int q; // down
	int sourceRate;
	int targetRate;

	int sourceDataLength;
	int targetDataLength;

	// continue==1 cache
	float *tailDataArr; // gcd
	int tailDataLength;

};

static void _resampleObj_calInterpArr(ResampleObj resampleObj);
static void _resampleObj_calInterpDeltaArr(ResampleObj resampleObj);

// NULL dataArr1
static float *_resampleObj_dealData(ResampleObj resampleObj,float *dataArr1,int dataLength1);
static void _resampleObj_resample(ResampleObj resampleObj,float *dataArr1,float *dataArr2);

/***
	qualType beat, use kaiser
	Best zeroNum 64 nbit 9 beta 14.7697 roll-off 0.9476
	Mid zeroNum 32 nbit 9 beta 11.6626  roll-off 0.8988
	Fast zeroNum 16 nbit 9 beta 8.5555	roll-off 0.85
****/
int resampleObj_new(ResampleObj *resampleObj,ResampleQualityType *qualType,int *isScale,int *isContinue){
	int status=0;

	ResampleQualityType _qualType=ResampleQuality_Best;

	int zeroNum=64; // 64 ->50
	int nbit=9; // 9 ->13

	WindowType winType=Window_Kaiser;
	float value=14.7696565; // 14.7696565 ->12.9846

	float rollOff=0.9475937; // 0.9475937 ->0.917347
	
	if(qualType){
		_qualType=*qualType;
	}

	if(_qualType==ResampleQuality_Mid){
		zeroNum=32;
		value=11.6625806;
		rollOff=0.8987969;
		nbit=9;
	}
	else if(_qualType==ResampleQuality_Fast){
		zeroNum=16;
		value=8.5555046;
		rollOff=0.85;
		nbit=9;
	}

	status=resampleObj_newWithWindow(resampleObj,
									&zeroNum,&nbit,
									&winType,&value,
									&rollOff,
									isScale,
									isContinue);

	return status;
}

/***
	sinc right
	zeroNum 64, 16/32/64
	nbit 9, 5~9 each zero-cross 1<<nbit samples
	winType hann
	value kaiser/gauss 5/2.5
	rollOff 0.945, 0.8~0.95
****/
int resampleObj_newWithWindow(ResampleObj *resampleObj,
							int *zeroNum,int *nbit,
							WindowType *winType,float *value,
							float *rollOff,
							int *isScale,
							int *isContinue){
	int status=0;

	int _isContinue=0;
	int _isScale=0;

	int _zeroNum=64;
	int _nbit=9;

	WindowType _winType=Window_Hann;
	float _value=0; // kaiser/gauss
	float _rollOff=0.945;

	ResampleObj resObj=NULL;

	int bitLength=0;
	int interpLength=0;

	if(isContinue){
		_isContinue=*isContinue;
	}

	if(isScale){
		_isScale=*isScale;
	}

	if(zeroNum){
		if(*zeroNum>0){
			_zeroNum=*zeroNum;
		}
	}

	if(nbit){
		if(*nbit>0&&*nbit<30){
			_nbit=*nbit;
		}
	}

	if(winType){
		if(*winType>Window_Rect){
			_winType=*winType;
		}
	}

	if(value){
		if(*value>=0){
			_value=*value;
		}

		if(_value==0){
			if(_winType==Window_Kaiser){
				_value=5;
			}
			else if(_winType==Window_Gauss){
				_value=2.5;
			}
		}
	}

	if(rollOff){
		if(*rollOff>0&&*rollOff<=1){
			_rollOff=*rollOff;
		}
	}
	
	resObj=*resampleObj=(ResampleObj )calloc(1, sizeof(struct OpaqueResample ));

	// 1. init
	bitLength=(1<<_nbit);
	interpLength=_zeroNum*bitLength+1;

	resObj->isContinue=_isContinue;
	resObj->isScale=_isScale;

	resObj->winType=_winType;
	resObj->value=_value;
	resObj->rollOff=_rollOff;

	resObj->zeroNum=_zeroNum;
	resObj->bitLength=bitLength;
	resObj->interpLength=interpLength;

	resObj->ratio=0.5;
	resObj->p=1;
	resObj->q=2;

	resObj->tailDataArr=__vnew(resObj->q+1, NULL);

	resObj->sourceRate=32000;
	resObj->targetRate=16000;

	// 2. cal interpArr
	_resampleObj_calInterpArr(resObj);

	// 3. cal interpDeltaArr
	resObj->interpDeltaArr=__vnew(interpLength, NULL);
	_resampleObj_calInterpDeltaArr(resObj);

	return status;
}

/***
	continue 0 
		sourceDataLength=dataLength
	conttinue 1
		sourceDataLength=dataLength-dataLength%q
****/
int resampleObj_calDataLength(ResampleObj resampleObj,int dataLength){
	int isContinue=0;

	float ratio=0;
	int p=0;
	int q=0;

	int sourceDataLength=0;
	int targetDataLength=0;

	isContinue=resampleObj->isContinue;

	ratio=resampleObj->ratio;
	p=resampleObj->p; // up
	q=resampleObj->q; // down

	if(!isContinue){
		sourceDataLength=dataLength;
		targetDataLength=floorf(dataLength*ratio);
	}
	else{
		if(q>1){ // down
			sourceDataLength=dataLength-dataLength%q;
			targetDataLength=sourceDataLength*p/q;
		}
	}

	resampleObj->sourceDataLength=sourceDataLength;
	resampleObj->targetDataLength=targetDataLength;

	return targetDataLength;

}

void resampleObj_setSamplate(ResampleObj resampleObj,int sourceRate,int targetRate){
	float ratio=0;
	int p=0;
	int q=0;

	float *interpArr=NULL; 
	int interpLength=0;

	int gcd=0;

	if(sourceRate==targetRate||
		sourceRate<=0||targetRate<=0){
		return;
	}

	interpArr=resampleObj->interpArr;
	interpLength=resampleObj->interpLength;

	gcd=util_gcd(sourceRate>targetRate?sourceRate:targetRate, sourceRate>targetRate?targetRate:sourceRate);
	p=targetRate/gcd; // up->target
	q=sourceRate/gcd; // down->source

	if(q>resampleObj->q){
		free(resampleObj->tailDataArr);
		resampleObj->tailDataArr=__vnew(q+1, NULL);
	}

	ratio=targetRate/(float )sourceRate;
	if(ratio!=resampleObj->ratio&&
		(resampleObj->ratio<1||ratio<1)){ // update interpArr/interpDeltaArr
		// 1. restore
		if(resampleObj->ratio<1){
			for(int i=0;i<interpLength;i++){
				interpArr[i]/=resampleObj->ratio;
			}
		}

		// 2. update ratio interpDeltaArr
		resampleObj->ratio=ratio;
		_resampleObj_calInterpDeltaArr(resampleObj);
	}

	resampleObj->ratio=ratio;
	resampleObj->sourceRate=sourceRate;
	resampleObj->targetRate=targetRate;

	resampleObj->p=p;
	resampleObj->q=q;
}

void resampleObj_setSamplateRatio(ResampleObj resampleObj,float ratio){
	float *interpArr=NULL; 
	int interpLength=0;
	
	if(ratio<0){
		return;
	}

	interpArr=resampleObj->interpArr;
	interpLength=resampleObj->interpLength;

	if(ratio!=resampleObj->ratio&&
		(resampleObj->ratio<1||ratio<1)){ // update interpArr/interpDeltaArr
		// 1. restore
		if(resampleObj->ratio<1){
			for(int i=0;i<interpLength;i++){
				interpArr[i]/=resampleObj->ratio;
			}
		}

		// 2. update ratio interpDeltaArr
		resampleObj->ratio=ratio;
		_resampleObj_calInterpDeltaArr(resampleObj);
	}

	resampleObj->ratio=ratio;

	resampleObj->p=0;
	resampleObj->q=0;
}

void resampleObj_enableContinue(ResampleObj resampleObj,int flag){

	if(!flag){ // not real-time change
		resampleObj->tailDataLength=0;
	}
	resampleObj->isContinue=flag;

}

/***
	1. dealData
	2. calDataLength
	3. resample
	4. tailData
	5. scale
****/
int resampleObj_resample(ResampleObj resampleObj,float *dataArr1,int dataLength1,float *dataArr2){
	int length=0;

	float *dealArr=NULL;
	float *curArr=NULL;

	int length1=0;
	int tailDataLength=0;

	// 1. dealData
	dealArr=_resampleObj_dealData(resampleObj,dataArr1,dataLength1);
	if(dealArr){ 
		curArr=dealArr;
		length1=dataLength1+resampleObj->tailDataLength;
	}
	else{
		curArr=dataArr1;
		length1=dataLength1;
	}

	// 2. calDataLength
	resampleObj_calDataLength(resampleObj,length1);

	// 3. resample
	_resampleObj_resample(resampleObj,curArr,dataArr2);

	// 4. tailData
	if(resampleObj->isContinue&&dealArr){ // real-time
		tailDataLength=length1-resampleObj->sourceDataLength;
		for(int i=0;i<tailDataLength;i++){
			resampleObj->tailDataArr[i]=dealArr[resampleObj->sourceDataLength+i];
		}

		resampleObj->tailDataLength=tailDataLength;
	}

	// 5. scale
	if(resampleObj->isScale){
		int _len=0;
		float _value=0;

		_len=resampleObj->targetDataLength;
		_value=sqrtf(resampleObj->ratio);
		for(int i=0;i<_len;i++){
			dataArr2[i]/=_value;
		}
	}

	length=resampleObj->targetDataLength;
	if(dealArr){
		free(dealArr);
	}
	return length;
}

// NULL ==> dataArr1
static float *_resampleObj_dealData(ResampleObj resampleObj,float *dataArr1,int dataLength1){
	float *dataArr2=NULL;

	int isContinue=0;
	float *tailDataArr=NULL; 
	int tailDataLength=0;

	isContinue=resampleObj->isContinue;
	tailDataArr=resampleObj->tailDataArr;
	tailDataLength=resampleObj->tailDataLength;
	if(isContinue&&tailDataLength){ 
		dataArr2=__vnew(dataLength1+tailDataLength, NULL);
		for(int i=0;i<tailDataLength;i++){
			dataArr2[i]=tailDataArr[i];
		}

		for(int i=0;i<dataLength1;i++){
			dataArr2[i+tailDataLength]=dataArr1[i];
		}
	}

	return dataArr2;
}

static void _resampleObj_resample(ResampleObj resampleObj,float *dataArr1,float *dataArr2){
	float ratio=0; // newSamplate/origSamplate

	int bitLength=0; // 1<<nbit
	float *interpArr=NULL; // ratio<1 *ration
	float *interpDeltaArr=NULL;
	int interpLength=0;
	
	int sourceDataLength=0;
	int targetDataLength=0;

	float scale=0;
	float t=0;
	int step=0;

	ratio=resampleObj->ratio;

	bitLength=resampleObj->bitLength;
	interpArr=resampleObj->interpArr;
	interpDeltaArr=resampleObj->interpDeltaArr;
	interpLength=resampleObj->interpLength;

	sourceDataLength=resampleObj->sourceDataLength;
	targetDataLength=resampleObj->targetDataLength;

	scale=(1.0>ratio?ratio:1.0);
	step=floorf(scale*bitLength);

	// {	
	// 	printf("interpArr is \n");
	// 	__vdebug(interpArr, interpLength, 1);
	// 	printf("\n");

	// 	printf("interpDeltaArr is \n");
	// 	__vdebug(interpDeltaArr, interpLength, 1);
	// 	printf("\n");
	// }

	for(int i=0;i<targetDataLength;i++){
		int leftLen=0;
		int rightLen=0;

		int n=0;
		float factor=0;
		float factorValue=0;
		int offset=0;

		float delta=0;
		float w=0;

		int _value1=0;
		int _value2=0;

		t=i*1.0/ratio;
		n=floorf(t);

		// 1. left cal
		factor=scale*(t-n);

		factorValue=factor*bitLength;
		offset=floorf(factorValue);
		delta=factorValue-offset;

		_value1=n+1;
		_value2=(interpLength-offset)/step;
		leftLen=(_value1>_value2?_value2:_value1);
		for(int j=0;j<leftLen;j++){
			w=interpArr[offset+j*step]+delta*interpDeltaArr[offset+j*step];
			dataArr2[i]+=w*dataArr1[n-j];
		}

		// 2. right cal
		// invert
		factor=scale-factor;

		factorValue=factor*bitLength;
		offset=floorf(factorValue);
		delta=factorValue-offset;

		_value1=sourceDataLength-n-1;
		_value2=(interpLength-offset)/step;
		rightLen=(_value1>_value2?_value2:_value1);
		for(int j=0;j<rightLen;j++){
			w=interpArr[offset+j*step]+delta*interpDeltaArr[offset+j*step];
			dataArr2[i]+=w*dataArr1[n+j+1];
		}

		// 3. update t
		// t+=1.0/ratio;

	}
}

static void _resampleObj_calInterpDeltaArr(ResampleObj resampleObj){
	float *interpArr=NULL;
	float *interpDeltaArr=NULL;

	int interpLength=0;
	float ratio=0;

	interpArr=resampleObj->interpArr;
	interpDeltaArr=resampleObj->interpDeltaArr;

	interpLength=resampleObj->interpLength;
	ratio=resampleObj->ratio;

	if(ratio<1){
		for(int i=0;i<interpLength;i++){
			interpArr[i]*=ratio;
		}
	}

	__vdiff(interpArr, interpLength, NULL,interpDeltaArr);
	interpDeltaArr[interpLength-1]=0; // diff ???
}

static void _resampleObj_calInterpArr(ResampleObj resampleObj){
	WindowType type=Window_Hann;
	float value=0;
	float rollOff=1;;

	int zeroNum=0;
	int interpLength=0; // zero*bitLength+1

	float *interpArr=NULL; // ratio<1 *ration

	float *winArr=NULL;
	int order=0;

	type=resampleObj->winType;
	value=resampleObj->value;
	rollOff=resampleObj->rollOff;

	zeroNum=resampleObj->zeroNum;
	interpLength=resampleObj->interpLength;

	if(type<=Window_Rect){
		type=Window_Hann;
	}
	
	// 1. sinc
	interpArr=__vlinspace(0, zeroNum, interpLength, 0);
	// == __vsinc_low
	for(int i=0;i<interpLength;i++){
		interpArr[i]*=rollOff;
	}
	__vsinc(interpArr, interpLength, interpArr);
	for(int i=0;i<interpLength;i++){
		interpArr[i]*=rollOff;
	}

	// 2. window
	order=(interpLength-1)*2;
	if(type==Window_Rect){
		float _value=1;

		winArr=__vnew(order+1, &_value);
	}
	else if(type==Window_Hann){
		winArr=window_createHann(order+1,0);
	}
	else if(type==Window_Hamm){
		winArr=window_createHamm(order+1,0);
	}
	else if(type==Window_Blackman){
		winArr=window_createBlackman(order+1,0);
	}
	else if(type==Window_Kaiser){
		winArr=window_createKaiser(order+1,0,&value);
	}
	else if(type==Window_Bartlett){
		winArr=window_createBartlett(order+1,0);
	}
	else if(type==Window_Triang){
		winArr=window_createTriang(order+1,0);
	}
	else if(type==Window_Flattop){
		winArr=window_createFlattop(order+1,0);
	}
	else if(type==Window_Gauss){
		winArr=window_createGauss(order+1,0,&value);
	}
	else if(type==Window_Blackman_Harris){
		winArr=window_createBlackmanHarris(order+1,0);
	}
	else if(type==Window_Blackman_Nuttall){
		winArr=window_createBlackmanNuttall(order+1,0);
	}
	else if(type==Window_Bartlett_Hann){
		winArr=window_createBartlettHann(order+1,0);
	}
	else if(type==Window_Bohman){
		winArr=window_createBohman(order+1,0);
	}
	else if(type==Window_Tukey){
		winArr=window_createTukey(order+1,0,&value);
	}

	// 3. sinc*window
	__vmul(interpArr, winArr+(interpLength-1), interpLength, interpArr);

	resampleObj->interpArr=interpArr;

	free(winArr);
}

void resampleObj_free(ResampleObj resampleObj){
	float *interpArr=NULL; // ratio<1 *ration
	float *interpDeltaArr=NULL;
	float *tailDataArr=NULL; // gcd

	if(!resampleObj){
		return;
	}

	interpArr=resampleObj->interpArr;
	interpDeltaArr=resampleObj->interpDeltaArr;
	tailDataArr=resampleObj->tailDataArr;

	free(interpArr);
	free(interpDeltaArr);
	free(tailDataArr);

	free(resampleObj);
}

void resampleObj_debug(ResampleObj resampleObj){

}














