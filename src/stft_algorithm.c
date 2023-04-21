// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"

#include "dsp/flux_window.h"
#include "dsp/fft_algorithm.h"

#include "stft_algorithm.h"

#ifdef HAVE_OMP
#include <omp.h>
#endif

// cpu core number
int __kernelNum=0;

typedef enum{
	STFTExec_None=0,

	STFTExec_STFT,
	STFTExec_ISTFT,

} STFTExecType;

struct OpaqueSTFT{
	STFTExecType execType;
	int useFlag;

	FFTObj fftObj;
	FFTObj *fftObjArr;

	WindowType windowType;
	float *windowDataArr;
	float *addDataArr;

	int fftLength; // y=fftLength ???
	int slideLength;
	int isContinue;

	int isPad; 
	PaddingPositionType positionType; // center
	PaddingModeType modeType; // constant zero

	float padValue1; // center/left/right constant
	float padValue2; // center constant

	float *tailDataArr;
	int tailDataLength;

	float *curDataArr;
	int curDataLength; // >=dataLength

	int timeLength; // x=timeLength

	float *realArr; // fft r cache

	// inverse 
	int methodType; // 0 'overlap-add' 1 'weight overlap-add'

	float *winArr1; // fftLength
	float *winArr2;
	
	int timeLength1;
	float *normArr; // (timeLength1-1)*slideLength+fftLength

	// int _weakFlag; // curDataArr weak flag
};

static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int isPad,int *timeLength,int *tailLength);

// status 0 error
static int __stftObj_dealData(STFTObj stftObj,float *dataArr,int dataLength);

static int __stftObj_dealPadData(STFTObj stftObj,float *dataArr,int dataLength,int tLen,float *curDataArr);

static void __stftObj_stft(STFTObj stftObj,float *dataArr,float *mRealArr,float *mImageArr);

static int __isCOA(float *winArr,int winLength,int overlapLength);

int stftObj_new(STFTObj *stftObj,int radix2Exp,WindowType *windowType,int *slideLength,int *isContinue){
	int status=0;

	int _isContinue=0;
	WindowType _windowType=Window_Rect;
	float *windowDataArr=NULL;
	float *addDataArr=NULL;

	STFTObj stft=NULL;
	FFTObj fftObj=NULL;

	if (__kernelNum==0){
        #ifdef HAVE_OMP
        __kernelNum=omp_get_max_threads()/2;
        if (__kernelNum==0){
            __kernelNum=1;
        }
        #else
        __kernelNum=1;
        #endif
	}

	FFTObj *fftObjArr=(FFTObj *)calloc(__kernelNum, sizeof(FFTObj ));

	int fftLength=0;
	int _slideLength=0;
	float *realArr=NULL;

	float *tailDataArr=NULL;

	if(radix2Exp<1||radix2Exp>30){
		status=-100;
		return status;
	}

	fftLength=1<<radix2Exp;
	if(windowType){
		_windowType=*windowType;
	}

	_slideLength=fftLength/4;
	if(slideLength){
		if(*slideLength>0){ // &&*slideLength<=fftLength ->support not overlap
			_slideLength=*slideLength;
		}
	}

	if(isContinue){
		_isContinue=*isContinue;
	}

	stft=*stftObj=(STFTObj )calloc(1, sizeof(struct OpaqueSTFT ));

	windowDataArr=window_calFFTWindow(_windowType,fftLength);
	addDataArr=__vnew(fftLength, NULL);
	fftObj_new(&fftObj, radix2Exp);

    #pragma omp parallel for
	for(int i=0;i<__kernelNum;i++){
		fftObj_new(fftObjArr+i, radix2Exp);
	}

	realArr=(float *)calloc(fftLength, sizeof(float ));
	tailDataArr=(float *)calloc(fftLength, sizeof(float ));

	stft->tailDataArr=tailDataArr;

	stft->windowType=_windowType;
	stft->windowDataArr=windowDataArr;
	stft->addDataArr=addDataArr;

	stft->fftObj=fftObj;
	stft->fftLength=fftLength;
	stft->slideLength=_slideLength;

	stft->realArr=realArr;
	stft->isContinue=_isContinue;

	stft->isPad=0;
	stft->positionType=PaddingPosition_Center;
	stft->modeType=PaddingMode_Constant;
	stft->fftObjArr=fftObjArr;

	return status;
}

// default 512 fftLength/4
void stftObj_setSlideLength(STFTObj stftObj,int slideLength){
	int fftLength=0;

	fftLength=stftObj->fftLength;
	if(slideLength>0){ // &&slideLength<=fftLength support not overlap
		stftObj->slideLength=slideLength;
	}
}

void stftObj_enableContinue(STFTObj stftObj,int flag){

	stftObj->isContinue=flag;
}

// center zero padding
void stftObj_enablePadding(STFTObj stftObj,int flag){

	stftObj->isPad=flag;
}

// value1 => constant/left
void stftObj_setPadding(STFTObj stftObj,
					PaddingPositionType *positionType,PaddingModeType *modeType,
					float *value1,float *value2){

	if(stftObj->isPad){
		if(positionType){
			stftObj->positionType=*positionType;
		}

		if(modeType){
			stftObj->modeType=*modeType;
		}

		if(value1){
			stftObj->padValue1=*value1;
		}

		if(value2){
			stftObj->padValue2=*value2;
		}
	}
}

void stftObj_useWindowDataArr(STFTObj stftObj,float *winDataArr){
	stftObj->useFlag=1;
	memcpy(stftObj->windowDataArr, winDataArr, sizeof(float )*stftObj->fftLength);
}

float *stftObj_getWindowDataArr(STFTObj stftObj){

	return stftObj->windowDataArr;
}

int stftObj_calTimeLength(STFTObj stftObj,int dataLength){
	int fftLength=0; 
	int slideLength=0;
	int tailDataLength=0;

	int isContinue=0;
	int isPad=0;

	int timeLength=0;

	fftLength=stftObj->fftLength;
	slideLength=stftObj->slideLength;
	tailDataLength=stftObj->tailDataLength;

	isContinue=stftObj->isContinue;
	isPad=stftObj->isPad;

	if(!isPad){
		if(isContinue){
			dataLength+=tailDataLength;
		}

		if(dataLength<fftLength){
			return 0;
		}

		timeLength=(dataLength-fftLength)/slideLength+1;
	}
	else{
		if(dataLength<=0){
			return 0;
		}

		timeLength=dataLength/slideLength+1;
	}

	return timeLength;
}

void stftObj_stft(STFTObj stftObj,float *dataArr,int dataLength,float *mRealArr,float *mImageArr){
	int status=0;

	if(!dataArr||dataLength<=0){
		return;
	}

	// 1. 处理相关数据
	if (stftObj->isPad || stftObj->isContinue) {
	    status=__stftObj_dealData(stftObj,dataArr,dataLength);
        if(!status){
            return;
        }
	}
	else{
	    stftObj->timeLength = stftObj_calTimeLength(stftObj, dataLength);
	}

	// 2. stft
	__stftObj_stft(stftObj,dataArr,mRealArr,mImageArr);
	
	stftObj->execType=STFTExec_STFT;

}

int stftObj_calDataLength(STFTObj stftObj,int timeLength){
	int fftLength=0; // y=fftLength ???
	int slideLength=0;

	int dataLength=0;

	fftLength=stftObj->fftLength;
	slideLength=stftObj->slideLength;

	dataLength=(timeLength-1)*slideLength+fftLength;

	return dataLength;
}

// type 0(deault) 'weight overlap-add' 1 'overlap-add'
void stftObj_istft(STFTObj stftObj,float *mRealArr,float *mImageArr,int timeLength1,int methodType,float *dataArr){
	FFTObj fftObj=NULL;
	float *windowDataArr=NULL;

	int fftLength=0; // y=fftLength ???
	int slideLength=0;

	int dataLength=0;

	float *winArr1=NULL; // fftLength
	float *winArr2=NULL;

	float *normArr=NULL; // (timeLength1-1)*slideLength+fftLength;

	float *realArr=NULL;
	float *imageArr=NULL;

	fftObj=stftObj->fftObj;
	windowDataArr=stftObj->windowDataArr;

	fftLength=stftObj->fftLength;
	slideLength=stftObj->slideLength;

	dataLength=(timeLength1-1)*slideLength+fftLength;

	normArr=stftObj->normArr;

	realArr=__vnew(fftLength, NULL);
	imageArr=__vnew(fftLength, NULL);

	// 1. winArr1/winArr2
	if(!stftObj->winArr1){
		float e=0;

		winArr1=__vnew(fftLength, NULL);
		winArr2=__vnew(fftLength, NULL);
		
		if(methodType==0){ // 'weight'
			e=1;
		}
		else{ // 'overlap-add' 
			e=0;
		}

		__vpow(windowDataArr, e, fftLength, winArr1);
		__vpow(windowDataArr, e+1, fftLength, winArr2);

		stftObj->winArr1=winArr1;
		stftObj->winArr2=winArr2;

		stftObj->methodType=methodType;
	}
	else{
		winArr1=stftObj->winArr1;
		winArr2=stftObj->winArr2;
		if(stftObj->methodType!=methodType){
			float e=0;

			if(methodType==0){ // 'weight'
				e=1;
			}
			else{ // 'overlap-add'
				e=0;
			}

			__vpow(windowDataArr, e, fftLength, winArr1);
			__vpow(windowDataArr, e+1, fftLength, winArr2);

			stftObj->methodType=methodType;
		}
	}

	// 2. normArr
	if(timeLength1>stftObj->timeLength1||
		stftObj->timeLength1>2*timeLength1){

		free(normArr);
		normArr=__vnew(dataLength, NULL);

		stftObj->normArr=normArr;
		stftObj->timeLength1=timeLength1;
	}

	// 3. ifft&overlap add
	memset(normArr, 0, sizeof(float )*dataLength);
	for(int i=0;i<timeLength1;i++){
		fftObj_ifft(fftObj, mRealArr+i*fftLength, mImageArr+i*fftLength, realArr, imageArr);
		
		for(int j=i*slideLength,k=0;j<i*slideLength+fftLength;j++,k++){
			dataArr[j]=dataArr[j]+realArr[k]*winArr1[k];
			normArr[j]=normArr[j]+winArr2[k];
		}
	}

	// 4. norm
	for(int i=0;i<dataLength;i++){
		if(normArr[i]<1e-6){
			normArr[i]=1;
		}
	}

	__vdiv(dataArr, normArr, dataLength, NULL);

	free(realArr);
	free(imageArr);
}

void stftObj_free(STFTObj stftObj){
	FFTObj fftObj=NULL;
	FFTObj *fftObjArr=NULL;

	float *windowDataArr=NULL;
	float *addDataArr=NULL;
	float *tailDataArr=NULL;
	float *curDataArr=NULL;

	float *realArr=NULL;

	float *winArr1=NULL; // fftLength
	float *winArr2=NULL;

	float *normArr=NULL;

	if(!stftObj){
		return;
	}

	fftObj=stftObj->fftObj;
	fftObjArr=stftObj->fftObjArr;

	windowDataArr=stftObj->windowDataArr;
	addDataArr=stftObj->addDataArr;
	tailDataArr=stftObj->tailDataArr;
	curDataArr=stftObj->curDataArr;

	realArr=stftObj->realArr;

	winArr1=stftObj->winArr1;
	winArr2=stftObj->winArr2;

	normArr=stftObj->normArr;

	fftObj_free(fftObj);

    #pragma omp parallel for
	for(int i=0;i<__kernelNum;i++){
		fftObj_free(fftObjArr[i]);
	}
	free(fftObjArr);

	free(windowDataArr);
	free(addDataArr);
	free(tailDataArr);
	free(curDataArr);

	free(realArr);

	free(winArr1);
	free(winArr2);

	free(normArr);

	free(stftObj);
}

/***
	1. 根据isContinue处理curDataArr
	2. 计算stftLength
	3. 处理tailDataArr
****/
static int __stftObj_dealData(STFTObj stftObj,float *dataArr,int dataLength){
	int status=1;

	int fftLength=0; 
	int slideLength=0;

	int isContinue=0;
	int isPad=0;

	float *tailDataArr=NULL;
	int tailDataLength=0;

	float *curDataArr=NULL;
	int curDataLength=0; 

	int timeLength=0;

	int timeLen=0;
	int tailLen=0;

	int totalLength=0;

	fftLength=stftObj->fftLength;
	slideLength=stftObj->slideLength;

	isContinue=stftObj->isContinue;
	isPad=stftObj->isPad;

	tailDataArr=stftObj->tailDataArr;
	tailDataLength=stftObj->tailDataLength;

	curDataArr=stftObj->curDataArr;
	curDataLength=stftObj->curDataLength;

	timeLength=stftObj->timeLength;

	// 1. curDataArr相关
	if(!isPad){ // 非填充
		if(isContinue){
			totalLength=tailDataLength+dataLength;
		}
		else{
			totalLength=dataLength;
		}

		if(totalLength<fftLength){
			tailLen=totalLength;
			status=0;
		}

		if(status){
			__calTimeAndTailLen(totalLength, fftLength, slideLength,isPad, &timeLen, &tailLen);
		}
	}
	else{ // 填充
		totalLength=dataLength;
		__calTimeAndTailLen(totalLength, fftLength, slideLength,isPad, &timeLen, &tailLen);
	}

	if(status){ // 存在timeLen 计算curDataArr相关
		if(totalLength>curDataLength||
			curDataLength>2*totalLength){

			free(curDataArr);
			curDataArr=(float *)calloc(totalLength+fftLength, sizeof(float ));
		}

		curDataLength=0;
		if(!isPad){ // 非填充fftLength
			if(isContinue&&tailDataLength<0){
				memcpy(curDataArr, dataArr-tailDataLength, (dataLength+tailDataLength)*sizeof(float ));
				curDataLength=(dataLength+tailDataLength);
			}
			else{
				if(isContinue&&tailDataLength>0){ // 连续且存在尾数据
					memcpy(curDataArr, tailDataArr, tailDataLength*sizeof(float ));
					curDataLength+=tailDataLength;
				}

				memcpy(curDataArr+curDataLength, dataArr, dataLength*sizeof(float ));
				curDataLength+=dataLength;
			}
		}
		else{ // 填充
			// tailLen=0; // padding is use tailLen
			curDataLength=__stftObj_dealPadData(stftObj,dataArr,dataLength,tailLen,curDataArr);
		}

		// 2. tailDataArr相关
		tailDataLength=0; // reset !!!
		if(isContinue&&!isPad){ // 连续非填充且存在tail
			if(tailLen>0){
				memcpy(tailDataArr,curDataArr+(curDataLength-tailLen),tailLen*sizeof(float ));
			}
			
			tailDataLength=tailLen;
		}
	}
	else{
		// 2. tailDataArr相关
		if(isContinue&&!isPad){ // 连续非填充且存在tail
			if(tailLen>0){
				if(tailDataLength>=0){
					memcpy(tailDataArr+tailDataLength,dataArr,dataLength*sizeof(float ));
				}
				else{
					memcpy(tailDataArr,dataArr-tailDataLength,(dataLength+tailDataLength)*sizeof(float ));
				}
			}
			
			tailDataLength=tailLen;
		}
		else{
			tailDataLength=0;
		}
	}

	// 3. 更新
	stftObj->tailDataLength=tailDataLength;

	stftObj->curDataArr=curDataArr;
	stftObj->curDataLength=curDataLength;

	stftObj->timeLength=timeLen;
	return status;
}

static int __stftObj_dealPadData(STFTObj stftObj,float *dataArr,int dataLength,int tLen,float *curDataArr){
	int fftLength=0; 
	int slideLength=0;

	int isContinue=0;
	int isPad=0;

	float *tailDataArr=NULL;
	int tailDataLength=0;

	int curDataLength=0; 

	float value1=0;
	float value2=0;

	int totalLength=0;
	int startIndex=0;

	int leftLen=0;
	int rightLen=0;

	fftLength=stftObj->fftLength;
	slideLength=stftObj->slideLength;

	isContinue=stftObj->isContinue;
	isPad=stftObj->isPad;

	tailDataArr=stftObj->tailDataArr;
	tailDataLength=stftObj->tailDataLength;

	value1=stftObj->padValue1;
	value2=stftObj->padValue2;

	// 1. 拼接数据
	if(stftObj->positionType==PaddingPosition_Center){
		startIndex=fftLength/2;
	}
	else if(stftObj->positionType==PaddingPosition_Left){
		startIndex=fftLength;
	}
	else if(stftObj->positionType==PaddingPosition_Right){
		startIndex=0;
	}

	// if(isContinue&&tailDataLength){ // 连续且存在尾数据
	// 	memcpy(curDataArr+startIndex, tailDataArr, tailDataLength*sizeof(float ));
	// 	curDataLength+=tailDataLength;
	// }
	
	if(dataLength>tLen){ // ??? 注意临界值问题
		memcpy(curDataArr+(startIndex+curDataLength), dataArr, (dataLength-tLen)*sizeof(float ));
	}
	curDataLength+=(dataLength-tLen);

	// 2. 填充数据
	leftLen=fftLength/2;
	rightLen=fftLength-fftLength/2;
	if(stftObj->modeType==PaddingMode_Constant){
		if(stftObj->positionType==PaddingPosition_Center){ // constant 两边填充
			__vpad_center1(curDataArr,curDataLength,leftLen,rightLen,value1,value2);
		}
		else if(stftObj->positionType==PaddingPosition_Left){ // 头部填充
			__vpad_left1(curDataArr,curDataLength,fftLength,value1);
		}
		else if(stftObj->positionType==PaddingPosition_Right){ // 尾部填充
			__vpad_right1(curDataArr,curDataLength,fftLength,value1);
		}
	}
	else if(stftObj->modeType==PaddingMode_Reflect&&curDataLength>1){ // reflect 反射
		if(stftObj->positionType==PaddingPosition_Center){ // 两边填充
			__vpad_center2(curDataArr,curDataLength,leftLen,rightLen);
		}
		else if(stftObj->positionType==PaddingPosition_Left){ // 开始填充
			__vpad_left2(curDataArr,curDataLength,fftLength);
		}
		else if(stftObj->positionType==PaddingPosition_Right){ // 尾部填充
			__vpad_right2(curDataArr,curDataLength,fftLength);
		}
	}
	else if(stftObj->modeType==PaddingMode_Wrap&&curDataLength>1){ // wrap 
		if(stftObj->positionType==PaddingPosition_Center){ // 两边填充
			__vpad_center3(curDataArr,curDataLength,leftLen,rightLen);
		}
		else if(stftObj->positionType==PaddingPosition_Left){ // 开始填充
			__vpad_left3(curDataArr,curDataLength,fftLength);
		}
		else if(stftObj->positionType==PaddingPosition_Right){ // 尾部填充
			__vpad_right3(curDataArr,curDataLength,fftLength);
		}
	}

	totalLength=curDataLength+fftLength;
	return totalLength;
}

static void __fft(STFTObj stftObj,FFTObj fftObj,int step,float *dataArr,float *mRealArr,float *mImageArr){
	int fftLength=0;
	int slideLength=0;

	float *windowDataArr=NULL;
	float *addDataArr=NULL;

	fftLength=stftObj->fftLength;
	slideLength=stftObj->slideLength;
	windowDataArr=stftObj->windowDataArr;
	addDataArr=__vnew(fftLength, NULL);

	for(int i=0;i<step;i++){
		__vmul(dataArr+i*slideLength, windowDataArr, fftLength, addDataArr);
//		__vdebug(addDataArr, fftLength, 1);
		fftObj_fft(fftObj,addDataArr,NULL,mRealArr+i*fftLength,mImageArr+i*fftLength);
	}

	free(addDataArr);
}

static void __stftObj_stft(STFTObj stftObj,float *dataArr, float *mRealArr,float *mImageArr){
	FFTObj fftObj=NULL;
	FFTObj *fftObjArr=NULL;
	int useFlag=0;

	WindowType windowType=Window_Rect;

	float *curDataArr=NULL;
	float *windowDataArr=NULL;
	float *addDataArr=NULL;

	int fftLength=0;
	int slideLength=0;
	int timeLength=0;// x=curSTFTCount

	int k=0;
	int block=0;
	int mod=0;

//	float *realArr=NULL;
    float *_arr = NULL;

	fftObj=stftObj->fftObj;
	fftObjArr=stftObj->fftObjArr;

	useFlag=stftObj->useFlag;

	windowType=stftObj->windowType;

	curDataArr=stftObj->curDataArr;
	windowDataArr=stftObj->windowDataArr;
	addDataArr=stftObj->addDataArr;

	fftLength=stftObj->fftLength;
	slideLength=stftObj->slideLength;
	timeLength=stftObj->timeLength;

//	realArr=stftObj->realArr;

	k=__kernelNum;
    block=timeLength/k;
    mod=timeLength%k;

    if (stftObj->isPad||stftObj->isContinue){
        _arr=curDataArr;
    }
    else{
        _arr=dataArr;
    }

    #ifdef HAVE_OMP
    if(timeLength==1){
        __fft(stftObj,fftObj,1,_arr,mRealArr,mImageArr);
    }else if(timeLength<__kernelNum){
        omp_set_num_threads(timeLength);

        #pragma omp parallel for
        for(int i=0;i<timeLength;i++){
            __fft(stftObj,fftObjArr[i],1,_arr+i*slideLength,mRealArr+i*fftLength,mImageArr+i*fftLength);
        }
    }else{
        omp_set_num_threads(__kernelNum);

        #pragma omp parallel for
        for(int i=0;i<k;i++){
            __fft(stftObj,fftObjArr[i],block,_arr+i*block*slideLength,mRealArr+i*block*fftLength,mImageArr+i*block*fftLength);
        }

        if(mod){
            __fft(stftObj,fftObj,mod,_arr+k*block*slideLength,mRealArr+k*block*fftLength,mImageArr+k*block*fftLength);
        }
    }
    #else
    for(int i=0;i<timeLength;i++){
        // memcpy(realArr, curDataArr+i*slideLength, sizeof(float )*fftLength);
        // fftObj_fft(fftObj,realArr,NULL,mRealArr+i*fftLength,mImageArr+i*fftLength);

        if(useFlag||windowType!=Window_Rect){
            __vmul(_arr+i*slideLength, windowDataArr, fftLength, addDataArr);
            fftObj_fft(fftObj,addDataArr,NULL,mRealArr+i*fftLength,mImageArr+i*fftLength);
        }
        else{
            fftObj_fft(fftObj,_arr+i*slideLength,NULL,mRealArr+i*fftLength,mImageArr+i*fftLength);
        }
    }
    #endif
}

/***
	isPad =0
		dataLength>=fftLength
		t=(dataLength-fftLength)/slideLength+1 
	isPad =1
		dataLength>0
		t=(dataLength+fftLength-fftLength)/slideLength+1 ;padding fftLength
****/
static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int isPad,int *timeLength,int *tailLength){
	int timeLen=0;
	int tailLen=0;

	if(!isPad){ // 非填充
		timeLen=(dataLength-fftLength)/slideLength+1;
		tailLen=(dataLength-fftLength)%slideLength+(fftLength-slideLength);
	}
	else{
		timeLen=dataLength/slideLength+1;
		if(timeLen>1){ // padding is use tailLen
			tailLen=dataLength%slideLength;
		}
	}

	if(timeLength){
		*timeLength=timeLen;
	}

	if(tailLength){
		*tailLength=tailLen;
	}
}

void stftObj_debug(STFTObj stftObj){
	int fftLength=0;
	int slideLength=0;

	int timeLength=0; // x=timeLength

	fftLength=stftObj->fftLength;
	slideLength=stftObj->slideLength;

	timeLength=stftObj->timeLength;

	printf("stft params is: fftLength=%d, slideLength=%d, timeLength=%d\n",fftLength,slideLength,timeLength);
}

static int __isCOA(float *winArr,int winLength,int overlapLength){
	int flag=0;


	return flag;
}















