// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"

#include "_pitch_yin.h"

struct OpaquePitchYIN{
	int isContinue;

	FFTObj fftObj;

	int fftLength;
	int slideLength;
	int autoLength; // autocorr length

	int minIndex; // min/maxFre
	int maxIndex;

	// yinLength <diffLength <fftLength
	int diffLength; // fftLength-autoLength
	int yinLength; // maxIndex-minIndex+1 <diffLength-1 
	int timeLength;

	// yin result ->timeLength
	float *mTroughArr; // timeLength*(yinLength/2+1)
	float *mFreArr;
	int *lenArr;

	float *mYinArr; // timeLength*yinLength
	float *mInterpArr; 

	float *mNumArr; //  timeLength*yinLength
	float *mDenArr;

	float *mDiffArr; // timeLength*diffLength
	float *mMeanArr; // timeLength*maxIndex

	// cache data ->fftLength
	float *realArr1; 
	float *imageArr1;

	float *realArr2; 
	float *imageArr2;

	float *realArr3; 
	float *imageArr3;

	float *energyArr1;
	float *energyArr2; 
	
	float *dataArr1;

	// continue相关数据
	float *tailDataArr;
	int tailDataLength;

	float *curDataArr;
	int curDataLength;

	int samplate;
	float thresh; // 0.1 good

	int isDebug;
};

static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int *timeLength,int *tailLength);

static void __pitchYINObj_initData(PitchYINObj pitchYINObj,int length);

static int __pitchYINObj_dealData(PitchYINObj pitchYINObj,float *dataArr,int dataLength);
static void __pitchYINObj_pitch(PitchYINObj pitchYINObj,float *freArr,float *troughArr,float *minArr);

static void __pitchYINObj_calDiff(PitchYINObj pitchYINObj);
static void __pitchYINObj_calInterp(PitchYINObj pitchYINObj);
static void __pitchYINObj_dealResult(PitchYINObj pitchYINObj,float *freArr,float *troughArr,float *minArr);

int pitchYINObj_new(PitchYINObj *pitchYINObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,int *autoLength,
				int *isContinue){
	int status=0;

	int _samplate=32000;
	float _lowFre=27; // 27.5
	float _highFre=2094; // 2093
	int _radix2Exp=12;
	int _slideLength=0;
	int _autoLength=0;
	int _isContinue=0;

	int fftLength=0;
	int minIndex=0; // min/maxFre
	int maxIndex=0;

	int diffLength=0;
	int yinLength=0;

	float thresh=0.1;

	FFTObj fftObj=NULL;
	PitchYINObj pitch=NULL;

	pitch=*pitchYINObj=(PitchYINObj )calloc(1,sizeof(struct OpaquePitchYIN ));

	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
		}
	}

	if(lowFre){
		if(*lowFre>=27){
			_lowFre=*lowFre;
		}
	}

	if(highFre){
		if(*highFre>_lowFre&&*highFre<_samplate/2){
			_highFre=*highFre;
		}
		else{
			_lowFre=27;
			_highFre=2093;
		}
	}

	if(radix2Exp){
		if(*radix2Exp>=1&&*radix2Exp<=30){
			_radix2Exp=*radix2Exp;
		}
	}

	fftLength=1<<_radix2Exp;
	_slideLength=fftLength/4;
	if(slideLength){
		if(*slideLength>0){ // &&*slideLength<=fftLength support not overlap
			_slideLength=*slideLength;
		}
	}

	_autoLength=fftLength/2;
	if(autoLength){
		if(*autoLength>=0&&*autoLength<fftLength){
			_autoLength=*autoLength;
		}
	}

	diffLength=fftLength-_autoLength;

	if(isContinue){
		_isContinue=*isContinue;
	}

	minIndex=floorf(_samplate/_highFre);
	maxIndex=ceilf(_samplate/_lowFre);
	if(maxIndex>diffLength-1){
		maxIndex=diffLength-1;
	}

	yinLength=maxIndex-minIndex+1;

	fftObj_new(&fftObj, _radix2Exp);
	
	pitch->isContinue=_isContinue;

	pitch->fftObj=fftObj;

	pitch->fftLength=fftLength;
	pitch->slideLength=_slideLength;
	pitch->autoLength=_autoLength;

	pitch->minIndex=minIndex;
	pitch->maxIndex=maxIndex;

	pitch->diffLength=diffLength;
	pitch->yinLength=yinLength;

	pitch->samplate=_samplate;
	pitch->thresh=thresh;

	// cache data ->fftLength
	__pitchYINObj_initData(pitch,fftLength);
	
	return status;
}

static void __pitchYINObj_initData(PitchYINObj pitchYINObj,int length){

	pitchYINObj->realArr1=__vnew(length, NULL);
	pitchYINObj->imageArr1=__vnew(length, NULL);

	pitchYINObj->realArr2=__vnew(length, NULL);
	pitchYINObj->imageArr2=__vnew(length, NULL);

	pitchYINObj->realArr3=__vnew(length, NULL);
	pitchYINObj->imageArr3=__vnew(length, NULL);

	pitchYINObj->energyArr1=__vnew(length, NULL);
	pitchYINObj->energyArr2=__vnew(length, NULL);

	pitchYINObj->dataArr1=__vnew(length, NULL);

	pitchYINObj->tailDataArr=__vnew(length, NULL);
}

// default 0.1 thresh>0&&thresh<1
void pitchYINObj_setThresh(PitchYINObj pitchYINObj,float thresh){

	// if(thresh>0&&thresh<1){
	// 	pitchYINObj->thresh=thresh;
	// }

	if(thresh>0){
		pitchYINObj->thresh=thresh;
	}
}

void pitchYINObj_pitch(PitchYINObj pitchYINObj,float *dataArr,int dataLength,
					float *freArr,float *valueArr1,float *valueArr2){
	int status=0;

	if(!dataArr||dataLength<=0){
		return;
	}

	// 1. 处理相关数据
	status=__pitchYINObj_dealData(pitchYINObj,dataArr,dataLength);
	if(!status){
		return;
	}

	// 2. yin
	__pitchYINObj_pitch(pitchYINObj,freArr,valueArr1,valueArr2);
}

int pitchYINObj_getTroughData(PitchYINObj pitchYINObj,float **mFreArr,float **mTroughArr,int **lenArr){
	int mLen=0;

	mLen=pitchYINObj->yinLength/2+1;
	if(mFreArr){
		*mFreArr=pitchYINObj->mFreArr;
	}

	if(mTroughArr){
		*mTroughArr=pitchYINObj->mTroughArr;
	}

	if(lenArr){
		*lenArr=pitchYINObj->lenArr;
	}

	return mLen;
}


static void __pitchYINObj_pitch(PitchYINObj pitchYINObj,float *freArr,float *troughArr,float *minArr){

	__pitchYINObj_calDiff(pitchYINObj);
	__pitchYINObj_calInterp(pitchYINObj);
	__pitchYINObj_dealResult(pitchYINObj,freArr,troughArr,minArr);
}

static void __pitchYINObj_calDiff(PitchYINObj pitchYINObj){
	FFTObj fftObj=NULL;

	int fftLength=0;
	int slideLength=0;
	int autoLength=0; // autocorr length

	int minIndex=0; // min/maxFre
	int maxIndex=0;

	int diffLength=0;
	int yinLength=0;
	int timeLength=0;

	float *mYinArr=NULL;
	float *mInterpArr=NULL; 

	float *mNumArr=NULL;
	float *mDenArr=NULL;

	float *mDiffArr=NULL;
	float *mMeanArr=NULL;

	float *realArr1=NULL; 
	float *imageArr1=NULL;

	float *realArr2=NULL; 
	float *imageArr2=NULL;

	float *realArr3=NULL; 
	float *imageArr3=NULL;

	float *energyArr1=NULL;
	float *energyArr2=NULL; 
	
	float *dataArr1=NULL; 

	float *curDataArr=NULL;

	fftObj=pitchYINObj->fftObj;

	fftLength=pitchYINObj->fftLength;
	slideLength=pitchYINObj->slideLength;
	autoLength=pitchYINObj->autoLength;

	minIndex=pitchYINObj->minIndex;
	maxIndex=pitchYINObj->maxIndex;

	diffLength=pitchYINObj->diffLength;
	yinLength=pitchYINObj->yinLength;
	timeLength=pitchYINObj->timeLength;

	mYinArr=pitchYINObj->mYinArr;
	mInterpArr=pitchYINObj->mInterpArr;

	mNumArr=pitchYINObj->mNumArr;
	mDenArr=pitchYINObj->mDenArr;

	mDiffArr=pitchYINObj->mDiffArr;
	mMeanArr=pitchYINObj->mMeanArr;

	realArr1=pitchYINObj->realArr1;
	imageArr1=pitchYINObj->imageArr1;

	realArr2=pitchYINObj->realArr2;
	imageArr2=pitchYINObj->imageArr2;

	realArr3=pitchYINObj->realArr3;
	imageArr3=pitchYINObj->imageArr3;

	energyArr1=pitchYINObj->energyArr1;
	energyArr2=pitchYINObj->energyArr2;

	dataArr1=pitchYINObj->dataArr1;

	curDataArr=pitchYINObj->curDataArr;

	for(int i=0;i<timeLength;i++){
		// 0. reset
		memset(realArr1, 0, sizeof(float )*fftLength);
		memset(imageArr1, 0, sizeof(float )*fftLength);

		memset(realArr2, 0, sizeof(float )*fftLength);
		memset(imageArr2, 0, sizeof(float )*fftLength);

		memset(realArr3, 0, sizeof(float )*fftLength);
		memset(imageArr3, 0, sizeof(float )*fftLength);

		// 1. auto correlation --> realArr1
		fftObj_fft(fftObj,curDataArr+i*slideLength,NULL,realArr1,imageArr1);

		for(int j=0;j<=autoLength;j++){
			dataArr1[j]=curDataArr[i*slideLength+autoLength-j];
		}
		fftObj_fft(fftObj,dataArr1,NULL,realArr2,imageArr2);

		__vcmul(realArr1, imageArr1, realArr2, imageArr2, fftLength, realArr3, imageArr3);

		memset(realArr1, 0, sizeof(float )*fftLength);
		memset(imageArr1, 0, sizeof(float )*fftLength);
		fftObj_ifft(fftObj, realArr3, imageArr3, realArr1, imageArr1);

		for(int j=autoLength,k=0;j<fftLength;j++,k++){
			float _value=0;

			_value=realArr1[j];
			if(fabs(_value)>=1e-6){
				realArr1[k]=_value;
			}
			else {
				realArr1[k]=0;
			}
		}

		// 2. energy -->energyArr2
		for(int j=0;j<fftLength;j++){
			float _value=0;

			_value=curDataArr[i*slideLength+j];
			if(j==0){
				energyArr1[j]=_value*_value;
			}
			else{
				energyArr1[j]=energyArr1[j-1]+_value*_value;
			}
		}

		for(int j=0;j<diffLength;j++){
			float _value=0;

			_value=energyArr1[autoLength+j]-energyArr1[j];
			if(fabs(_value)>=1e-6){
				energyArr2[j]=_value;
			}
			else {
				energyArr2[j]=0;
			}
		}

		// 3. difference
		for(int j=0;j<diffLength;j++){ // timeLength*diffLength
			mDiffArr[i*diffLength+j]=energyArr2[0]+energyArr2[j]-2*realArr1[j];
		}
	}

	// 4. cumu mean norm
	for(int i=0;i<timeLength;i++){ // mNumArr ->timeLength*yinLength
		for(int j=minIndex,k=0;j<maxIndex+1;j++,k++){
			mNumArr[i*(maxIndex-minIndex+1)+k]=mDiffArr[i*diffLength+j];
		}
	}

	for(int i=0;i<timeLength;i++){ // mMeanArr ->timeLength*maxIndex
		for(int j=1,k=0;j<maxIndex+1;j++,k++){
			float _value=0;

			_value=mDiffArr[i*diffLength+j];
			if(k==0){
				mMeanArr[i*maxIndex+k]=_value;
			}
			else{
				mMeanArr[i*maxIndex+k]=mMeanArr[i*maxIndex+k-1]+_value;
			}
		}

		for(int j=1,k=0;j<maxIndex+1;j++,k++){
			mMeanArr[i*maxIndex+k]/=j;
		}
	}

	for(int i=0;i<timeLength;i++){ // mDenArr ->timeLength*yinLength
		for(int j=minIndex-1,k=0;j<maxIndex;j++,k++){
			mDenArr[i*(maxIndex-minIndex+1)+k]=mMeanArr[i*maxIndex+j];
		}
	}

	for(int i=0;i<timeLength;i++){
		for(int j=0;j<yinLength;j++){
			mYinArr[i*yinLength+j]=mNumArr[i*yinLength+j]/(mDenArr[i*yinLength+j]+1e-16);
		}
	}

	// {
	// 	printf("mMeanArr is:\n");
	// 	__mdebug(mMeanArr, timeLength, maxIndex, 1);
	// 	printf("\n");
	// }
}

static void __pitchYINObj_calInterp(PitchYINObj pitchYINObj){
	int yinLength=0; // maxLength-minLength+1
	int timeLength=0;

	float *mYinArr=NULL;
	float *mInterpArr=NULL;

	float num=0;
	float den=0;

	float offset=0;

	float value1=0;
	float value2=0;
	float value3=0;

	yinLength=pitchYINObj->yinLength;
	timeLength=pitchYINObj->timeLength;

	mYinArr=pitchYINObj->mYinArr;
	mInterpArr=pitchYINObj->mInterpArr;

	memset(mInterpArr, 0, sizeof(float )*timeLength*yinLength);
	for(int i=0;i<timeLength;i++){
		for(int j=0;j<yinLength-2;j++){
			value1=mYinArr[i*yinLength+j];
			value2=mYinArr[i*yinLength+j+1];
			value3=mYinArr[i*yinLength+j+2];

			num=(value3-value1)/2;
			den=(value1+value3-2*value2)/2;
			
			offset=-num/(2*den+1e-16);
			if(fabsf(offset)<=1){ // 2*p ∈[-1,1]
				mInterpArr[i*yinLength+j+1]=offset;
			}
			else{
				mInterpArr[i*yinLength+j+1]=0;
			}
		}
	}
}

static void __pitchYINObj_dealResult(PitchYINObj pitchYINObj,float *freArr,float *troughArr,float *minArr){
	int yinLength=0; // maxLength-minLength+1
	int timeLength=0;

	float *mTroughArr=NULL; // timeLength*(yinLength/2+1)
	float *mFreArr=NULL;
	int *lenArr=NULL;

	float *mYinArr=NULL;
	float *mInterpArr=NULL;

	int minIndex=0;

	int samplate=0;
	float thresh=0;

	float *arr1=NULL;
	float offset=0;

	int troughIndex=0;

	yinLength=pitchYINObj->yinLength;
	timeLength=pitchYINObj->timeLength;

	mTroughArr=pitchYINObj->mTroughArr;
	mFreArr=pitchYINObj->mFreArr;
	lenArr=pitchYINObj->lenArr;

	mYinArr=pitchYINObj->mYinArr;
	mInterpArr=pitchYINObj->mInterpArr;

	minIndex=pitchYINObj->minIndex;

	samplate=pitchYINObj->samplate;
	thresh=pitchYINObj->thresh;
	
	for(int i=0;i<timeLength;i++){
		arr1=mYinArr+i*yinLength;
		
		// 1. trough < thresh
		troughIndex=-1;
		for(int j=0;j<yinLength-1;j++){
			if(j==0){
				if(arr1[j]<arr1[j+1]&&
					arr1[j]<thresh){

					troughIndex=j;
					if(troughArr){
						troughArr[i]=arr1[j];
					}

					break;
				}
			}
			else{
				if(arr1[j]<=arr1[j+1]&&
					arr1[j]<arr1[j-1]&&
					arr1[j]<thresh){

					troughIndex=j;
					if(troughArr){
						troughArr[i]=arr1[j];
					}

					break;
				}
			}
		}

		// 2. cal fre
		if(troughIndex!=-1){ 
			offset=mInterpArr[i*yinLength+troughIndex];
			freArr[i]=samplate/(minIndex+troughIndex+offset);
		}
		
		if(minArr){
			__vmin(arr1, yinLength, minArr+i);
		}
	}

	// cal mTrough/mFreArr/lenArr
	{
		int len=0;
		int mLen=0;

		mLen=yinLength/2+1;
		for(int i=0;i<timeLength;i++){
			len=0;
			arr1=mYinArr+i*yinLength;
			
			for(int j=0;j<yinLength-1;j++){
				if(j==0){
					if(arr1[j]<arr1[j+1]&&
						arr1[j]<thresh){

						mTroughArr[i*mLen+len]=arr1[j];

						offset=mInterpArr[i*yinLength+j];
						mFreArr[i*mLen+len]=samplate/(minIndex+j+offset);

						len++;
					}
				}
				else{
					if(arr1[j]<=arr1[j+1]&&
						arr1[j]<arr1[j-1]&&
						arr1[j]<thresh){
						
						mTroughArr[i*mLen+len]=arr1[j];

						offset=mInterpArr[i*yinLength+j];
						mFreArr[i*mLen+len]=samplate/(minIndex+j+offset);

						len++;
					}
				}
			}

			lenArr[i]=len;
		}
	}

	// trough vector and min vector
	if(pitchYINObj->isDebug){
		float *minArr=NULL;
		int *mIndexArr=NULL;

		float *troughArr=NULL;
		int *tIndexArr=NULL;
		float *freArr=NULL;

		int len=0;

		minArr=__vnew(timeLength, NULL);
		mIndexArr=__vnewi(timeLength, NULL);

		troughArr=__vnew(yinLength/2+1, NULL);
		tIndexArr=__vnewi(yinLength/2+1, NULL);
		freArr=__vnew(yinLength/2+1, NULL);

		printf("pitch debug start ------\n");

		for(int i=0;i<timeLength;i++){
			len=0;
			arr1=mYinArr+i*yinLength;
			
			for(int j=0;j<yinLength-1;j++){
				if(j==0){
					if(arr1[j]<arr1[j+1]&&
						arr1[j]<thresh){

						tIndexArr[len]=j;
						troughArr[len]=arr1[j];

						offset=mInterpArr[i*yinLength+j];
						freArr[len]=samplate/(minIndex+j+offset);

						len++;
					}
				}
				else{
					if(arr1[j]<=arr1[j+1]&&
						arr1[j]<arr1[j-1]&&
						arr1[j]<thresh){
						
						tIndexArr[len]=j;
						troughArr[len]=arr1[j];

						offset=mInterpArr[i*yinLength+j];
						freArr[len]=samplate/(minIndex+j+offset);

						len++;
					}
				}
			}

			mIndexArr[i]=__vmin(arr1, yinLength, minArr+i);

			printf("index[%d]:%.3f=\n",i,i*1.0*pitchYINObj->slideLength/pitchYINObj->samplate);

			if(len){
				printf(" 	");
				__vdebug(troughArr, len, 0);

				printf(" 	");
				__vdebugi(tIndexArr, len, 0);

				printf(" 	");
				__vdebug(freArr, len, 0);
			}

			printf(" 	");
			printf("min:%d ,%f \n",mIndexArr[i],minArr[i]);
		}

		printf("\n");
		printf("minArr:\n");

		printf(" 	");
		__vdebug(minArr, timeLength, 0);

		printf(" 	");
		__vdebugi(mIndexArr, timeLength, 0);

		printf("pitch debug end ------\n");

		free(minArr);
		free(mIndexArr);

		free(troughArr);
		free(tIndexArr);
		free(freArr);
	}
}

int pitchYINObj_calTimeLength(PitchYINObj pitchYINObj,int dataLength){
	int fftLength=0; 
	int slideLength=0;
	int tailDataLength=0;

	int isContinue=0;

	int timeLength=0;

	fftLength=pitchYINObj->fftLength;
	slideLength=pitchYINObj->slideLength;
	tailDataLength=pitchYINObj->tailDataLength;

	isContinue=pitchYINObj->isContinue;

	if(isContinue){
		dataLength+=tailDataLength; // outTimeLength
	}

	if(dataLength<fftLength){
		return 0;
	}

	timeLength=(dataLength-fftLength)/slideLength+1;
	return timeLength;
}

void pitchYINObj_enableDebug(PitchYINObj pitchYINObj,int isDebug){

	pitchYINObj->isDebug=isDebug;
}

void pitchYINObj_free(PitchYINObj pitchYINObj){

	if(pitchYINObj){
		fftObj_free(pitchYINObj->fftObj);

		free(pitchYINObj->mTroughArr);
		free(pitchYINObj->mFreArr);
		free(pitchYINObj->lenArr);

		free(pitchYINObj->mYinArr);
		free(pitchYINObj->mInterpArr);

		free(pitchYINObj->mNumArr);
		free(pitchYINObj->mDenArr);

		free(pitchYINObj->mDiffArr);
		free(pitchYINObj->mMeanArr);

		free(pitchYINObj->realArr1);
		free(pitchYINObj->imageArr1);

		free(pitchYINObj->realArr2);
		free(pitchYINObj->imageArr2);

		free(pitchYINObj->realArr3);
		free(pitchYINObj->imageArr3);

		free(pitchYINObj->energyArr1);
		free(pitchYINObj->energyArr2);

		free(pitchYINObj->dataArr1);

		free(pitchYINObj->tailDataArr);
		free(pitchYINObj->curDataArr);

		free(pitchYINObj);
	}
}

static int __pitchYINObj_dealData(PitchYINObj pitchYINObj,float *dataArr,int dataLength){
	int status=1;

	int fftLength=0; 
	int slideLength=0;

	int isContinue=0;

	float *tailDataArr=NULL;
	int tailDataLength=0;

	float *curDataArr=NULL;
	int curDataLength=0; 

	int maxIndex=0;

	int diffLength=0;
	int yinLength=0;
	int timeLength=0;

	int timeLen=0;
	int tailLen=0;

	int totalLength=0;

	fftLength=pitchYINObj->fftLength;
	slideLength=pitchYINObj->slideLength;

	isContinue=pitchYINObj->isContinue;

	tailDataArr=pitchYINObj->tailDataArr;
	tailDataLength=pitchYINObj->tailDataLength;

	curDataArr=pitchYINObj->curDataArr;
	curDataLength=pitchYINObj->curDataLength;

	maxIndex=pitchYINObj->maxIndex;

	diffLength=pitchYINObj->diffLength;
	yinLength=pitchYINObj->yinLength;
	timeLength=pitchYINObj->timeLength;

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
		__calTimeAndTailLen(totalLength, fftLength, slideLength, &timeLen, &tailLen);
	}

	if(status){ // 存在timeLen 计算curDataArr相关
		if(totalLength>curDataLength||
			curDataLength>2*totalLength){

			free(curDataArr);
			curDataArr=(float *)calloc(totalLength+fftLength, sizeof(float ));
		}

		curDataLength=0;
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

		// tailDataArr相关
		tailDataLength=0; // reset !!!
		if(isContinue){ 
			if(tailLen>0){
				memcpy(tailDataArr,curDataArr+(curDataLength-tailLen),tailLen*sizeof(float ));
			}
			
			tailDataLength=tailLen;
		}

		// yin 相关
		if(pitchYINObj->timeLength<timeLen||
			pitchYINObj->timeLength>timeLen*2){ // 更新缓存
			free(pitchYINObj->mTroughArr);
			free(pitchYINObj->mFreArr);
			free(pitchYINObj->lenArr);

			free(pitchYINObj->mYinArr);
			free(pitchYINObj->mInterpArr);

			free(pitchYINObj->mNumArr);
			free(pitchYINObj->mDenArr);

			free(pitchYINObj->mDiffArr);
			free(pitchYINObj->mMeanArr);

			pitchYINObj->mTroughArr=__vnew(timeLen*(yinLength/2+1),NULL);
			pitchYINObj->mFreArr=__vnew(timeLen*(yinLength/2+1),NULL);
			pitchYINObj->lenArr=__vnewi(timeLen,NULL);

			pitchYINObj->mYinArr=__vnew(timeLen*yinLength,NULL);
			pitchYINObj->mInterpArr=__vnew(timeLen*yinLength,NULL);

			pitchYINObj->mNumArr=__vnew(timeLen*yinLength,NULL);
			pitchYINObj->mDenArr=__vnew(timeLen*yinLength,NULL);

			pitchYINObj->mDiffArr=__vnew(timeLen*diffLength,NULL);
			pitchYINObj->mMeanArr=__vnew(timeLen*maxIndex,NULL);
		}
	}
	else{
		if(isContinue){ 
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

	pitchYINObj->tailDataLength=tailDataLength;

	pitchYINObj->curDataArr=curDataArr;
	pitchYINObj->curDataLength=curDataLength;

	pitchYINObj->timeLength=timeLen;

	return status;
}

static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int *timeLength,int *tailLength){
	int timeLen=0;
	int tailLen=0;

	timeLen=(dataLength-fftLength)/slideLength+1;
	tailLen=(dataLength-fftLength)%slideLength+(fftLength-slideLength);

	if(timeLength){
		*timeLength=timeLen;
	}

	if(tailLength){
		*tailLength=tailLen;
	}
}









