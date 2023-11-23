// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"

#include "_pitch_pef.h"

struct OpaquePitchPEF{
	int isContinue;

	FFTObj fftObj1; // 2*fftLength
	FFTObj fftObj2; // (4||8)*fftLength ???

	int fftLength; 
	int slideLength;
	int radix2Exp; // fftLength

	int xcorrFFTLength; // (4||8)*fftLength ???

	int timeLength;

	int minIndex; // edge
	int maxIndex;

	float *winDataArr; // fftLength

	// freband
	float *linearFreBandArr; // fftLength+1
	float *logFreBandArr; // fftLength*2
	float *bandWidthArr; // fftLength*2

	// estimate filter
	float alpha; // 10 
	float beta; // 0.5 0~1
	float gamma; // 1.8 >1

	float *filterArr; // fftLength*8 ->xcorrFFTLength
	int filterPadNum;

	float *mPowerArr; // timeLength*(2*fftLength)
	float *mInterpArr; // timeLength*(8*fftLength) ->xcorrFFTLength
	float *mXcorrArr; // timeLength*(8*fftLength) ->xcorrFFTLength

	// cache data
	float *realArr1; // fftLength*8 ->xcorrFFTLength
	float *imageArr1;

	float *realArr2; 
	float *imageArr2;

	float *realArr3;
	float *imageArr3;

	float *dataArr1; // fftLength*2

	// continue
	float *tailDataArr; // fftLength
	int tailDataLength;

	float *curDataArr;
	int curDataLength;

	int samplate;
	WindowType winType;

	float lowFre;
	float highFre;
	float cutFre; // >highFre

	int isDebug;
};

static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int *timeLength,int *tailLength);

static void __pitchPEFObj_initData(PitchPEFObj pitchPEFObj);
static void __pitchPEFObj_calEstimateFilter(PitchPEFObj pitchPEFObj);

static int __pitchPEFObj_dealData(PitchPEFObj pitchPEFObj,float *dataArr,int dataLength);

static void __pitchPEFObj_calInterp(PitchPEFObj pitchPEFObj);
static void __pitchPEFObj_calXcorr(PitchPEFObj pitchPEFObj);
static void __pitchPEFObj_dealResult(PitchPEFObj pitchPEFObj,float *freArr);

/***
	samplate 32000
	lowFre 32,
	highFre 2000
	cutFre 4000, >highFre

	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	alpha 10 >0, beta 0.5 0~1, gamma 1.8 >1
	
	isContinue 0
****/
int pitchPEFObj_new(PitchPEFObj *pitchPEFObj,
				int *samplate,float *lowFre,float *highFre,float *cutFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				float *alpha,float *beta,float *gamma,
				int *isContinue){
	int status=0;
	
	float _lowFre=32; // 27.5
	float _highFre=2000; // 2093/4186
	float _cutFre=4000;

	int _samplate=32000;
	int _radix2Exp=12;
	int _slideLength=0;
	WindowType _winType=Window_Hamm;
	int _isContinue=0;

	float _alpha=10; // 10 >1
	float _beta=0.5; // 0.5 0~1
	float _gamma=1.8; // 1.8 >1

	int fftLength=0;

	FFTObj fftObj1=NULL;
	
	PitchPEFObj pitch=NULL;

	pitch=*pitchPEFObj=(PitchPEFObj )calloc(1,sizeof(struct OpaquePitchPEF ));

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
			_lowFre=32;
			_highFre=2000;
		}
	}

	if(cutFre){
		if(*cutFre>=_highFre){
			_cutFre=*cutFre;
		}
		else{
			_cutFre=_highFre;
		}
	}

	if(radix2Exp){
		if(*radix2Exp>=1&&*radix2Exp<=30){
			_radix2Exp=*radix2Exp;
		}
	}

	if(windowType){
		_winType=*windowType;
	}

	if(alpha){
		if(*alpha>0){
			_alpha=*alpha;
		}
	}

	if(beta){
		if(*beta>0){
			_beta=*beta;
		}
	}

	if(gamma){
		if(*gamma>1){
			_gamma=*gamma;
		}
	}

	fftLength=1<<_radix2Exp;
	_slideLength=fftLength/4;
	if(slideLength){
		if(*slideLength>0){ // &&*slideLength<=fftLength support not overlap
			_slideLength=*slideLength;
		}
	}

	if(isContinue){
		_isContinue=*isContinue;
	}

	fftObj_new(&fftObj1, _radix2Exp+1);

	pitch->fftObj1=fftObj1;

	pitch->fftLength=fftLength;
	pitch->slideLength=_slideLength;
	pitch->radix2Exp=_radix2Exp;

	pitch->alpha=_alpha;
	pitch->beta=_beta;
	pitch->gamma=_gamma;

	pitch->samplate=_samplate;
	pitch->winType=_winType;

	pitch->lowFre=_lowFre;
	pitch->highFre=_highFre;
	pitch->cutFre=_cutFre;

	pitch->isContinue=_isContinue;

	__pitchPEFObj_initData(pitch);
	__pitchPEFObj_calEstimateFilter(pitch);

	return status;
}

void pitchPEFObj_pitch(PitchPEFObj pitchPEFObj,float *dataArr,int dataLength,
					float *freArr){
	int status=0;

	if(!dataArr||dataLength<=0){
		return;
	}

	// 1. deal data
	status=__pitchPEFObj_dealData(pitchPEFObj,dataArr,dataLength);
	if(!status){
		return;
	}

	// 2. power&interp
	__pitchPEFObj_calInterp(pitchPEFObj);

	// 3. cross correlation
	__pitchPEFObj_calXcorr(pitchPEFObj);

	// 4. peak pick
	__pitchPEFObj_dealResult(pitchPEFObj,freArr);

}

static void __pitchPEFObj_calInterp(PitchPEFObj pitchPEFObj){
	FFTObj fftObj1=NULL;

	int fftLength=0;
	int slideLength=0;
	int timeLength=0;

	float *winDataArr=NULL; // fftLength

	// freband
	float *linearFreBandArr=NULL; // fftLength+1
	float *logFreBandArr=NULL; // fftLength*2
	float *bandWidthArr=NULL; // fftLength*2

	int xcorrFFTLength=0;
	int filterPadNum=0;

	float *mPowerArr=NULL; 
	float *mInterpArr=NULL; 

	float *realArr1=NULL; 
	float *imageArr1=NULL;

	float *dataArr1=NULL; 
	float *curDataArr=NULL;

	fftObj1=pitchPEFObj->fftObj1;

	fftLength=pitchPEFObj->fftLength;
	slideLength=pitchPEFObj->slideLength;
	timeLength=pitchPEFObj->timeLength;

	winDataArr=pitchPEFObj->winDataArr;

	linearFreBandArr=pitchPEFObj->linearFreBandArr;
	logFreBandArr=pitchPEFObj->logFreBandArr;
	bandWidthArr=pitchPEFObj->bandWidthArr;

	xcorrFFTLength=pitchPEFObj->xcorrFFTLength;
	filterPadNum=pitchPEFObj->filterPadNum;

	mPowerArr=pitchPEFObj->mPowerArr;
	mInterpArr=pitchPEFObj->mInterpArr;

	realArr1=pitchPEFObj->realArr1;
	imageArr1=pitchPEFObj->imageArr1;

	dataArr1=pitchPEFObj->dataArr1;
	curDataArr=pitchPEFObj->curDataArr;

	for(int i=0;i<timeLength;i++){
		// 0. reset
		memset(dataArr1, 0, sizeof(float )*fftLength*2);

		// 1. fft
		memcpy(dataArr1, curDataArr+i*slideLength, sizeof(float )*fftLength);
		__vmul(dataArr1, winDataArr, fftLength, NULL);

		fftObj_fft(fftObj1, dataArr1, NULL, realArr1, imageArr1);
		for(int j=0;j<fftLength+1;j++){
			mPowerArr[i*fftLength*2+j]=realArr1[j]*realArr1[j]+imageArr1[j]*imageArr1[j];
		}

		// 2. interp
		__vinterp_linear(linearFreBandArr, mPowerArr+i*fftLength*2, fftLength+1, logFreBandArr, fftLength*2, mInterpArr+(i*xcorrFFTLength+filterPadNum));

		// 3. weight
		memset(mInterpArr+i*xcorrFFTLength, 0, sizeof(float )*filterPadNum);
		__vmul(mInterpArr+(i*xcorrFFTLength+filterPadNum), bandWidthArr, fftLength*2, NULL);
	}
}

static void __pitchPEFObj_calXcorr(PitchPEFObj pitchPEFObj){
	FFTObj fftObj2=NULL; 
	int timeLength=0;

	float *filterArr=NULL;
	int filterPadNum=0;

	int xcorrFFTLength=0; 

	float *mInterpArr=NULL;
	float *mXcorrArr=NULL; 

	float *realArr1=NULL;
	float *imageArr1=NULL;

	float *realArr2=NULL;
	float *imageArr2=NULL;

	float *realArr3=NULL;
	float *imageArr3=NULL;

	fftObj2=pitchPEFObj->fftObj2;
	timeLength=pitchPEFObj->timeLength;

	filterArr=pitchPEFObj->filterArr;
	filterPadNum=pitchPEFObj->filterPadNum;

	xcorrFFTLength=pitchPEFObj->xcorrFFTLength;

	mInterpArr=pitchPEFObj->mInterpArr;
	mXcorrArr=pitchPEFObj->mXcorrArr;

	realArr1=pitchPEFObj->realArr1;
	imageArr1=pitchPEFObj->imageArr1;

	realArr2=pitchPEFObj->realArr2;
	imageArr2=pitchPEFObj->imageArr2;

	realArr3=pitchPEFObj->realArr3;
	imageArr3=pitchPEFObj->imageArr3;

	fftObj_fft(fftObj2, filterArr, NULL, realArr2, imageArr2);
	for(int i=0;i<xcorrFFTLength;i++){
		imageArr2[i]=-imageArr2[i];
	}

	for(int i=0;i<timeLength;i++){
		fftObj_fft(fftObj2, mInterpArr+i*xcorrFFTLength, NULL, realArr3, imageArr3);
		__vcmul(realArr3, imageArr3, realArr2, imageArr2, xcorrFFTLength, realArr3, imageArr3);

		fftObj_ifft(fftObj2, realArr3, imageArr3, mXcorrArr+i*xcorrFFTLength, imageArr1);
	}
}

static void __pitchPEFObj_dealResult(PitchPEFObj pitchPEFObj,float *freArr){
	int timeLength=0;

	int fftLength=0;
	int xcorrFFTLength=0; // (4||8)*fftLength ???
	int filterPadNum=0;

	int minIndex=0; // edge
	int maxIndex=0;

	float *logFreBandArr=NULL; 
	float *mXcorrArr=NULL;

	float *realArr3=NULL;

	float value1=0;
	int index1=0;

	int len=0;

	timeLength=pitchPEFObj->timeLength;

	fftLength=pitchPEFObj->fftLength;
	xcorrFFTLength=pitchPEFObj->xcorrFFTLength;
	filterPadNum=pitchPEFObj->filterPadNum;

	minIndex=pitchPEFObj->minIndex;
	maxIndex=pitchPEFObj->maxIndex;

	logFreBandArr=pitchPEFObj->logFreBandArr;
	mXcorrArr=pitchPEFObj->mXcorrArr;

	realArr3=pitchPEFObj->realArr3;

	len=(maxIndex<fftLength*2+filterPadNum-1?maxIndex+1:fftLength*2+filterPadNum-1);
	for(int i=0;i<timeLength;i++){
		memcpy(realArr3,mXcorrArr+(i*xcorrFFTLength+(xcorrFFTLength-len)),sizeof(float )*len);
		memcpy(realArr3+len,mXcorrArr+i*xcorrFFTLength,sizeof(float )*(len+1));

		util_peakPick(realArr3+(maxIndex+1), 2*len-maxIndex, minIndex, maxIndex, 1, 1, &value1, &index1);
		freArr[i]=logFreBandArr[index1];
	}
}

static void __pitchPEFObj_initData(PitchPEFObj pitchPEFObj){
	int fftLength=0;

	int minIndex=-1; // edge
	int maxIndex=0;

	float *winDataArr=NULL; // fftLength

	// freband
	float *linearFreBandArr=NULL; // fftLength+1
	float *logFreBandArr=NULL; // fftLength*2
	float *bandWidthArr=NULL; 

	int samplate=0;
	WindowType winType=Window_Hamm;

	float lowFre=0;
	float highFre=0;
	float cutFre=0;

	float fre1=0;

	fftLength=pitchPEFObj->fftLength;

	samplate=pitchPEFObj->samplate;
	winType=pitchPEFObj->winType;

	lowFre=pitchPEFObj->lowFre;
	highFre=pitchPEFObj->highFre;
	cutFre=pitchPEFObj->cutFre;

	winDataArr=window_calFFTWindow(winType, fftLength);
	linearFreBandArr=__vlinspace(0, samplate/2, fftLength+1, 0);
	
	fre1=(samplate/2>cutFre?cutFre:samplate/2-1);
	logFreBandArr=__vlogspace(1, log10f(fre1), fftLength*2, 0);

	for(int i=1;i<fftLength*2;i++){
		if(highFre<logFreBandArr[i]){
			if(logFreBandArr[i]-highFre<highFre-logFreBandArr[i-1]){
				maxIndex=i;
			}
			else{
				maxIndex=i-1;
			}

			break;
		}

		if(minIndex!=-1){
			continue;
		}

		if(lowFre<logFreBandArr[i]){
			if(logFreBandArr[i]-lowFre<lowFre-logFreBandArr[i-1]){
				minIndex=i;
			}
			else{
				minIndex=i-1;
			}
		}
	}

	bandWidthArr=__vnew(fftLength*2, NULL);
	for(int i=2,j=1;i<fftLength*2;i++,j++){
		bandWidthArr[j]=(logFreBandArr[i]-logFreBandArr[i-2])/(2*fftLength*2);
	}

	bandWidthArr[0]=bandWidthArr[1];
	bandWidthArr[fftLength*2-1]=bandWidthArr[fftLength*2-2];

	pitchPEFObj->minIndex=minIndex;
	pitchPEFObj->maxIndex=maxIndex;

	pitchPEFObj->winDataArr=winDataArr;

	pitchPEFObj->linearFreBandArr=linearFreBandArr;
	pitchPEFObj->logFreBandArr=logFreBandArr;
	pitchPEFObj->bandWidthArr=bandWidthArr;

	pitchPEFObj->filterArr=__vnew(fftLength*8, NULL);

	pitchPEFObj->realArr1=__vnew(fftLength*8, NULL);
	pitchPEFObj->imageArr1=__vnew(fftLength*8, NULL);

	pitchPEFObj->realArr2=__vnew(fftLength*8, NULL);
	pitchPEFObj->imageArr2=__vnew(fftLength*8, NULL);

	pitchPEFObj->realArr3=__vnew(fftLength*8, NULL);
	pitchPEFObj->imageArr3=__vnew(fftLength*8, NULL);

	pitchPEFObj->dataArr1=__vnew(fftLength*2, NULL);

	pitchPEFObj->tailDataArr=__vnew(fftLength, NULL);
}

static int __pitchPEFObj_dealData(PitchPEFObj pitchPEFObj,float *dataArr,int dataLength){
	int status=1;

	int fftLength=0; 
	int slideLength=0;

	int isContinue=0;

	float *tailDataArr=NULL;
	int tailDataLength=0;

	float *curDataArr=NULL;
	int curDataLength=0; 

	int timeLength=0;

	int xcorrFFTLength=0; // (4||8)*fftLength ???

	float *mPowerArr=NULL; // timeLength*fftLength
	float *mInterpArr=NULL; // timeLength*(2*fftLength)
	float *mXcorrArr=NULL; // timeLength*(8*fftLength)

	int timeLen=0;
	int tailLen=0;

	int totalLength=0;

	fftLength=pitchPEFObj->fftLength;
	slideLength=pitchPEFObj->slideLength;

	isContinue=pitchPEFObj->isContinue;

	tailDataArr=pitchPEFObj->tailDataArr;
	tailDataLength=pitchPEFObj->tailDataLength;

	curDataArr=pitchPEFObj->curDataArr;
	curDataLength=pitchPEFObj->curDataLength;

	timeLength=pitchPEFObj->timeLength;
	xcorrFFTLength=pitchPEFObj->xcorrFFTLength;

	mPowerArr=pitchPEFObj->mPowerArr;
	mInterpArr=pitchPEFObj->mInterpArr;
	mXcorrArr=pitchPEFObj->mXcorrArr;

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

	if(status){ // has timeLen, cal curDataArr
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
			if(isContinue&&tailDataLength>0){ // has & tail
				memcpy(curDataArr, tailDataArr, tailDataLength*sizeof(float ));
				curDataLength+=tailDataLength;
			}

			memcpy(curDataArr+curDataLength, dataArr, dataLength*sizeof(float ));
			curDataLength+=dataLength;
		}

		// tailDataArr
		tailDataLength=0; // reset !!!
		if(isContinue){ 
			if(tailLen>0){
				memcpy(tailDataArr,curDataArr+(curDataLength-tailLen),tailLen*sizeof(float ));
			}
			
			tailDataLength=tailLen;
		}

		// update cache
		if(pitchPEFObj->timeLength<timeLen||
			pitchPEFObj->timeLength>timeLen*2){ 
			free(pitchPEFObj->mPowerArr);
			free(pitchPEFObj->mInterpArr);
			free(pitchPEFObj->mXcorrArr);
	
			pitchPEFObj->mPowerArr=__vnew(timeLen*fftLength*2,NULL);
			pitchPEFObj->mInterpArr=__vnew(timeLen*fftLength*8,NULL);
			pitchPEFObj->mXcorrArr=__vnew(timeLen*fftLength*8,NULL);
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

	pitchPEFObj->tailDataLength=tailDataLength;

	pitchPEFObj->curDataArr=curDataArr;
	pitchPEFObj->curDataLength=curDataLength;

	pitchPEFObj->timeLength=timeLen;

	return status;
}

int pitchPEFObj_calTimeLength(PitchPEFObj pitchPEFObj,int dataLength){
	int fftLength=0; 
	int slideLength=0;
	int tailDataLength=0;

	int isContinue=0;

	int timeLength=0;

	fftLength=pitchPEFObj->fftLength;
	slideLength=pitchPEFObj->slideLength;
	tailDataLength=pitchPEFObj->tailDataLength;

	isContinue=pitchPEFObj->isContinue;

	if(isContinue){
		dataLength+=tailDataLength; // outTimeLength
	}

	if(dataLength<fftLength){
		return 0;
	}

	timeLength=(dataLength-fftLength)/slideLength+1;
	return timeLength;
}

void pitchPEFObj_setFilterParams(PitchPEFObj pitchPEFObj,float alpha,float beta,float gamma){

	if(!(alpha>0&&beta>0&&gamma>1)){
		return;
	}

	if(alpha!=pitchPEFObj->alpha||beta!=pitchPEFObj->beta||gamma!=pitchPEFObj->gamma){
		__pitchPEFObj_calEstimateFilter(pitchPEFObj);
	}
}

static void __pitchPEFObj_calEstimateFilter(PitchPEFObj pitchPEFObj){
	int fftLength=0;

	FFTObj fftObj2=NULL; 
	int radix2Exp=0;

	int xcorrFFTLength=0;

	float *filterArr=NULL;
	int filterPadNum=0;

	float alpha=0; // 10 
	float beta=0; // 0.5 0~1
	float gamma=0; // 1.8 >1

	float *qArr=NULL;
	float *hArr=NULL;
	float *dArr=NULL;

	float value1=0;
	float value2=0;
	float det=0;

	fftLength=pitchPEFObj->fftLength;

	alpha=pitchPEFObj->alpha;
	beta=pitchPEFObj->beta;
	gamma=pitchPEFObj->gamma;

	filterArr=pitchPEFObj->filterArr;

	qArr=__vlogspace(log10f(beta), log10f(alpha+beta), fftLength, 0);
	hArr=__vnew(fftLength, NULL);
	dArr=__vnew(fftLength+1, NULL);

	// 1. hArr,filterPadNum
	for(int i=0;i<fftLength;i++){
		if(qArr[i]<1){
			filterPadNum++;
		}

		hArr[i]=1/(gamma-cosf(2*M_PI*qArr[i]));
	}

	// 2. diff
	dArr[0]=qArr[0];
	for(int i=1;i<fftLength;i++){
		dArr[i]=(qArr[i-1]+qArr[i])/2;
	}
	dArr[fftLength]=qArr[fftLength-1];

	for(int i=1;i<fftLength+1;i++){
		dArr[i-1]=dArr[i]-dArr[i-1];
	}

	// 3. sum
	value1=__vsum(dArr, fftLength);
	__vmul(dArr, hArr, fftLength, NULL);
	value2=__vsum(dArr, fftLength);
	det=value2/value1;

	// 4. xcorrFFTLength
	radix2Exp=pitchPEFObj->radix2Exp;
	if(filterPadNum){
		radix2Exp+=3;
	}
	else{
		radix2Exp+=2;
	}

	xcorrFFTLength=1<<radix2Exp;

	// 5. filterArr
	for(int i=0;i<fftLength;i++){
		filterArr[i]=hArr[i]-det;
	}

	// 6. fftObj
	fftObj_new(&fftObj2, radix2Exp);
	fftObj_free(pitchPEFObj->fftObj2);

	free(qArr);
	free(hArr);
	free(dArr);

	pitchPEFObj->fftObj2=fftObj2;
	pitchPEFObj->xcorrFFTLength=xcorrFFTLength;

	pitchPEFObj->filterPadNum=filterPadNum;
}

void pitchPEFObj_enableDebug(PitchPEFObj pitchPEFObj,int isDebug){

	pitchPEFObj->isDebug=isDebug;
}

void pitchPEFObj_free(PitchPEFObj pitchPEFObj){

	if(pitchPEFObj){
		fftObj_free(pitchPEFObj->fftObj1);
		fftObj_free(pitchPEFObj->fftObj2);

		free(pitchPEFObj->winDataArr);

		free(pitchPEFObj->linearFreBandArr);
		free(pitchPEFObj->logFreBandArr);
		free(pitchPEFObj->bandWidthArr);

		free(pitchPEFObj->filterArr);

		free(pitchPEFObj->mPowerArr);
		free(pitchPEFObj->mInterpArr);
		free(pitchPEFObj->mXcorrArr);

		free(pitchPEFObj->realArr1);
		free(pitchPEFObj->imageArr1);

		free(pitchPEFObj->realArr2);
		free(pitchPEFObj->imageArr2);

		free(pitchPEFObj->realArr3);
		free(pitchPEFObj->imageArr3);

		free(pitchPEFObj->dataArr1);

		free(pitchPEFObj->tailDataArr);
		free(pitchPEFObj->curDataArr);

		free(pitchPEFObj);
	}
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








