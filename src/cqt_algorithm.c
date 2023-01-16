// clang -g -c 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/fft_algorithm.h"
#include "dsp/dct_algorithm.h"
#include "dsp/resample_algorithm.h"

#include "filterbank/cqt_filterBank.h"
#include "filterbank/chroma_filterBank.h"

#include "stft_algorithm.h"
#include "cqt_algorithm.h"

struct OpaqueCQT{
	int isContinue; // 针对cqt 非stft/resample 
	int vFlag; // vqt标志

	STFTObj stftObj;
	ResampleObj resampleObj;

	int fftLength; // fftLength,timeLength,num
	int timeLength; // cqt timeLength isContinue 0 ==stft timeLength
	int stftTimeLength;

	// filterBank相关
	int num;
	int octaveNum;
	int binPerOctave;

	float *mRealFilterBankArr; // num*(fftLength/2+1) -> binPerOctave*(fftLength/2+1)
	float *mImageFilterBankArr;

	float *freBandArr; // num
	float *lenArr; // binPerOctave

	float *sLenArr; // 针对 num sqrt
	float *dLenArr; // 针对 octaveNum sqrt

	int samplate;
	WindowType windowType;
	SpectralFilterBankNormalType normType;
	int isScale;

	float minFre; // C1 32.703196
	int chromaNum; // 默认12
	float *mChromaFilterBank; // chromaNum*num 
	float *mSArr; // timeLength*num

	// continue相关数据
	float *tailDataArr;
	int tailDataLength;

	float *validDataArr;
	int validDataLength;

	// stft相关
	float *mRealArr1; // stft r,i(timeLength*fftLength) -> power r(timeLength*(fftLength/2+1)) 
	float *mImageArr1;

	float *mRealArr2; // timeLength*binPerOctave
	float *mImageArr2;

	int radix2Exp; // stft
	int slideLength;

	// dct相关
	FFTObj fftObj; // dct加速
	DCTObj dctObj; // 直接矩阵

	// deconv相关
	FFTObj devFFTObj;
	int devFFTLength;

	float *devDataArr; // spectral mag&fft mag

	float *devRealArr1; // fft
	float *devImageArr1;

	float *devRealArr2; // ifft
	float *devImageArr2;

	int isDebug;
};

static void _cqtObj_dealFilterBank(CQTObj cqtObj,int num,float minFre,int samplate,
								int binPerOctave,SpectralFilterBankNormalType normType,WindowType winType,
								float factor,float beta,float thresh,int vFlag);
static void _cqtObj_dealStft(CQTObj cqtObj,int fftLength,int slideLength,int isContinue);
static void _cqtObj_dealResample(CQTObj cqtObj);

static void _cqtObj_dealDeconv(CQTObj cqtObj);

static void _cqtObj_dealDCT(CQTObj cqtObj,int num);

// 处理tail/cur Data
static int _cqtObj_dealData(CQTObj cqtObj,float *dataArr,int dataLength);
static void _cqtObj_cqt(CQTObj cqtObj,float *dataArr,int dataLength,float *mRealArr,float *mImageArr);

static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int isContinue,int *timeLength,int *tailLength);

static int _calRadix2(int length);

int cqtObj_new(CQTObj *cqtObj,int num,int samplate,float minFre,int *isContinue){
	int status=0;

	status=cqtObj_newWith(cqtObj,num,
						&samplate,&minFre,NULL,
						NULL,NULL,NULL,
						NULL,NULL,isContinue,
						NULL,NULL);

	return status;
}

// 7*binPerOctave --> binPerOctave-fftLength 12-512/24-1024/36-1024/48-2048
int cqtObj_newWith(CQTObj *cqtObj,int num,
				int *samplate,float *minFre,int *binPerOctave,
				float *factor,float *beta,float *thresh,
				WindowType *windowType,int *slideLength,int *isContinue,
				SpectralFilterBankNormalType *normalType,int *isScale){
	int status=0;
	CQTObj cqt=NULL;

	int _samplate=32000;
	float _minFre=32.703196; // C1
	int _binPerOctave=12;

	float _factor=1;
	float _beta=0;
	float _thresh=0.01; // paper 0.0054
	
	int _slideLength=0;
	int _isContinue=0;
	int vFlag=0;

	float *tailDataArr=NULL;

	WindowType _windowType=Window_Hann;
	SpectralFilterBankNormalType _normalType=SpectralFilterBankNormal_None;

	int _isScale=1;

	cqt=*cqtObj=(CQTObj )calloc(1, sizeof(struct OpaqueCQT ));

	if(binPerOctave){
		if(*binPerOctave>0){
			_binPerOctave=*binPerOctave;
		}
	}

	if(_binPerOctave%12!=0){
		printf("binPerOctave is error\n");
		return -1;
	}

	if(num<_binPerOctave||num%_binPerOctave!=0){
		printf("num is error\n");
		return -1;
	}

	if(samplate){
		if(*samplate>0){
			_samplate=*samplate;
		}
	}

	if(minFre){
		if(*minFre>0){
			_minFre=*minFre;
		}
	}

	if(factor){
		if(*factor>0){
			_factor=*factor;
		}
	}

	if(beta){
		if(*beta>0){
			_beta=*beta;
		}
		
		if(_beta!=0){
			vFlag=1;
		}
	}

	if(thresh){
		if(*thresh>0){
			_thresh=*thresh;
		}
	}

	if(windowType){
		_windowType=*windowType;
	}

	if(slideLength){
		if(*slideLength>0){
			_slideLength=*slideLength;
		}
	}

	if(isContinue){
		_isContinue=*isContinue;
	}

	if(normalType){
		_normalType=*normalType;
	}

	if(isScale){
		_isScale=*isScale;
	}

	_cqtObj_dealFilterBank(cqt,num,_minFre,_samplate,
						_binPerOctave,_normalType,_windowType,
						_factor,_beta,_thresh,vFlag);

	if(_slideLength<=0){ // 支持非overlap >fftLength
		_slideLength=cqt->fftLength/4;
	}

	tailDataArr=(float *)calloc(cqt->fftLength, sizeof(float ));

	_cqtObj_dealStft(cqt,cqt->fftLength,_slideLength,_isContinue);
	_cqtObj_dealResample(cqt);
	_cqtObj_dealDCT(cqt,num);

	cqt->tailDataArr=tailDataArr;

	cqt->isScale=_isScale;
	cqt->isContinue=_isContinue;
	cqt->vFlag=vFlag;

	cqt->minFre=_minFre;

	return status;
}

static void _cqtObj_dealDCT(CQTObj cqtObj,int num){
	FFTObj fftObj=NULL;
	DCTObj dctObj=NULL;
	int r=0;

	r=_calRadix2(num);
	if(r){ // dct加速
		fftObj_new(&fftObj, r);
	}
	else{ // 直接dct
		dctObj_new(&dctObj, num,NULL);
	}

	cqtObj->fftObj=fftObj;
	cqtObj->dctObj=dctObj;
}

int cqtObj_calTimeLength(CQTObj cqtObj,int dataLength){
	int fftLength=0; 
	int slideLength=0;
	int tailDataLength=0;

	int isContinue=0;

	int timeLength=0;

	fftLength=cqtObj->fftLength;
	slideLength=cqtObj->slideLength;
	tailDataLength=cqtObj->tailDataLength;

	isContinue=cqtObj->isContinue;

	if(isContinue){
		dataLength+=tailDataLength; // outTimeLength

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

/***
	isContinue =1
		dataLength>=fftLength
		t=(dataLength-fftLength)/slideLength+1 
	isContinue =0
		dataLength>0
		t=(dataLength+fftLength-fftLength)/slideLength+1 ;padding fftLength
****/
static void __calTimeAndTailLen(int dataLength,int fftLength,int slideLength,int isContinue,int *timeLength,int *tailLength){
	int timeLen=0;
	int tailLen=0;

	if(isContinue){ // 连续
		timeLen=(dataLength-fftLength)/slideLength+1;
		tailLen=(dataLength-fftLength)%slideLength+(fftLength-slideLength);
	}
	else{
		timeLen=dataLength/slideLength+1;
		// if(timeLen>1){
		// 	tailLen=dataLength%slideLength;
		// }
	}

	if(timeLength){
		*timeLength=timeLen;
	}

	if(tailLength){
		*tailLength=tailLen;
	}
}

int cqtObj_getFFTLength(CQTObj cqtObj){
	int length=0;

	length=cqtObj->fftLength;
	return length;
}

float *cqtObj_getFreBandArr(CQTObj cqtObj){

	return cqtObj->freBandArr;
}

// isContinue 1 时处理tail/cur Data 
static int _cqtObj_dealData(CQTObj cqtObj,float *dataArr,int dataLength){
	int status=1;

	int isContinue=0;

	int fftLength=0; // fftLength,timeLength,num
	int timeLength=0;

	int slideLength=0;

	float *tailDataArr=NULL;
	int tailDataLength=0;

	float *validDataArr=NULL;
	int validDataLength=0;

	int timeLen=0;
	int tailLen=0;

	int totalLength=0;

	isContinue=cqtObj->isContinue;

	fftLength=cqtObj->fftLength;
	timeLength=cqtObj->timeLength;

	slideLength=cqtObj->slideLength;

	tailDataArr=cqtObj->tailDataArr;
	tailDataLength=cqtObj->tailDataLength;

	validDataArr=cqtObj->validDataArr;
	validDataLength=cqtObj->validDataLength;

	if(isContinue){ 
		totalLength=tailDataLength+dataLength;

		if(totalLength<fftLength){
			tailLen=totalLength;
			status=0;
		}

		if(status){
			__calTimeAndTailLen(totalLength, fftLength, slideLength,isContinue, &timeLen, &tailLen);
		}
	}
	else{
		totalLength=dataLength;
		__calTimeAndTailLen(totalLength, fftLength, slideLength,isContinue, &timeLen, &tailLen);
	}

	if(status){
		if(totalLength>validDataLength||
			validDataLength>2*totalLength){

			free(validDataArr);
			validDataArr=(float *)calloc(totalLength+fftLength, sizeof(float ));
		}

		validDataLength=0;
		if(isContinue&&tailDataLength<0){
			memcpy(validDataArr, dataArr-tailDataLength, (dataLength+tailDataLength)*sizeof(float ));
			validDataLength=(dataLength+tailDataLength);
		}
		else{
			if(isContinue&&tailDataLength>0){ // 连续且存在尾数据
				memcpy(validDataArr, tailDataArr, tailDataLength*sizeof(float ));
				validDataLength+=tailDataLength;
			}

			memcpy(validDataArr+validDataLength, dataArr, dataLength*sizeof(float ));
			validDataLength+=dataLength;
		}
		
		// 2. tailDataArr相关
		tailDataLength=0;
		if(isContinue){
			if(tailLen>0){
				memcpy(tailDataArr, validDataArr+(validDataLength-tailLen), tailLen*sizeof(float ));
			}

			tailDataLength=tailLen;
		}
	}
	else{
		// 2. tailDataArr相关
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

	cqtObj->tailDataLength=tailDataLength;

	cqtObj->validDataArr=validDataArr;
	cqtObj->validDataLength=validDataLength;

	cqtObj->timeLength=timeLen;
	return status;
}

void cqtObj_setScale(CQTObj cqtObj,int flag){

	cqtObj->isScale=flag;
}

void cqtObj_cqt(CQTObj cqtObj,float *dataArr,int dataLength,float *mRealArr3,float *mImageArr3){
	int status=0;

	if(!dataArr||dataLength<=0){
		return;
	}

	// 1. 处理相关数据
	status=_cqtObj_dealData(cqtObj,dataArr,dataLength);
	if(!status){
		return;
	}
	
	// 2. cqt
	_cqtObj_cqt(cqtObj,cqtObj->validDataArr,cqtObj->validDataLength,mRealArr3,mImageArr3);
}

/***
	cqt->chroma
	num(84) --> binPerOctave(12) --> chromaNum
****/
void cqtObj_chroma(CQTObj cqtObj,int *chromaNum,SpectralDataType *dataType,ChromaDataNormalType *normType,
				float *mRealArr1,float *mImageArr1,
				float *mDataArr3){
	int _chromaNum=12;

	SpectralDataType _dataType=SpectralData_Power;
	ChromaDataNormalType _normType=ChromaDataNormal_Max;

	int timeLength=0;

	int num=0;
	int binPerOctave=0;

	float *mChromaFilterBank=NULL;
	float *mSArr=NULL;

	float p1=1;
	int type1=0;

	if(chromaNum){
		_chromaNum=*chromaNum;
	}

	if(dataType){
		// if(*dataType<SpectralData_DB){
		// 	_dataType=*dataType;
		// }
		_dataType=*dataType;
	}

	if(normType){
		_normType=*normType;
	}

	timeLength=cqtObj->timeLength;

	num=cqtObj->num;
	binPerOctave=cqtObj->binPerOctave;

	if(_chromaNum>binPerOctave||binPerOctave%_chromaNum!=0){
		printf("chromaNum and binPerOctave not map!!!");
		return;
	}

	mChromaFilterBank=cqtObj->mChromaFilterBank;
	mSArr=cqtObj->mSArr;
	if(_chromaNum!=cqtObj->chromaNum){
		free(mSArr);
		free(mChromaFilterBank);
		
		mSArr=__vnew(timeLength*num, NULL);
		mChromaFilterBank=__vnew(_chromaNum*num, NULL);
		
		chroma_cqtFilterBank(_chromaNum,num,binPerOctave,
							&cqtObj->minFre,
							mChromaFilterBank);
	}

	for(int i=0;i<timeLength;i++){
		for(int j=0;j<num;j++){
			float _value=0;

			_value=mRealArr1[i*num+j]*mRealArr1[i*num+j]+mImageArr1[i*num+j]*mImageArr1[i*num+j];
			if(_dataType==SpectralData_Mag){
				_value=sqrtf(_value);
			}

			mSArr[i*num+j]=_value;
		}
	}

	__mdot1(mSArr,mChromaFilterBank,
		timeLength,num,
		_chromaNum,num,
		mDataArr3);

	// for(int i=0;i<timeLength;i++){
	// 	for(int j=0;j<num;j++){
	// 		float _value=0;

	// 		_value=mRealArr1[i*num+j]*mRealArr1[i*num+j]+mImageArr1[i*num+j]*mImageArr1[i*num+j];
	// 		if(_dataType==SpectralData_Mag){
	// 			_value=sqrtf(_value);
	// 		}

	// 		// sum == chromaMatrix.dot(S)
	// 		mDataArr3[i*binPerOctave+j%binPerOctave]+=_value;
	// 		// if(mDataArr3[i*binPerOctave+j%binPerOctave]<_value){ // max
	// 		// 	mDataArr3[i*binPerOctave+j%binPerOctave]=_value;
	// 		// }
	// 	}
	// }

	if(_normType!=ChromaDataNormal_None){
		if(_normType==ChromaDataNormal_Max){
			type1=1;
		}
		else if(_normType==ChromaDataNormal_Min){
			type1=2;
		}
		else{
			type1=0;
			if(_normType==ChromaDataNormal_P2){
				p1=2;
			}
		}

		__mnormalize(mDataArr3,timeLength,_chromaNum,1,type1,p1,mDataArr3);
	}

	cqtObj->chromaNum=_chromaNum;
	cqtObj->mChromaFilterBank=mChromaFilterBank;
	cqtObj->mSArr=mSArr;
}

void cqtObj_cqcc(CQTObj cqtObj,float *mDataArr1,int mLength,CepstralRectifyType *rectifyType,float *mDataArr2){
	FFTObj fftObj=NULL;
	DCTObj dctObj=NULL;

	CepstralRectifyType recType=CepstralRectify_Log;

	int num=0;
	int timeLength=0;

	float *mRealArr=NULL; 
	float *mImageArr=NULL;

	num=cqtObj->num;
	timeLength=cqtObj->timeLength;

	mRealArr=cqtObj->mRealArr1;
	mImageArr=cqtObj->mImageArr1;

	fftObj=cqtObj->fftObj;
	dctObj=cqtObj->dctObj;

	if(mLength>num){
		return;
	}

	if(rectifyType){
		recType=*rectifyType;
	}

	if(recType==CepstralRectify_CubicRoot){
		for(int i=0;i<timeLength*num;i++){
			mRealArr[i]=powf(mDataArr1[i], 1.0/3);
		}
	}
	else { // log
		for(int i=0;i<timeLength*num;i++){
			float _value=0;

			_value=mDataArr1[i];
			if(_value<1e-8){
				_value=1e-8;
			}

			mRealArr[i]=log10f(_value); // matlab canonic
		}
	}

	for(int i=0;i<timeLength;i++){
		if(fftObj){
			fftObj_dct(fftObj, mRealArr+i*num, mImageArr+i*num, 1);
		}
		else{
			dctObj_dct(dctObj, mRealArr+i*num,1,mImageArr+i*num);
		}
	}

	for(int i=0;i<timeLength;i++){
		for(int j=0;j<mLength;j++){
			mDataArr2[i*mLength+j]=mImageArr[i*num+j];
		}
	}
}

void cqtObj_cqhc(CQTObj cqtObj,float *mDataArr1,int hcNum,float *mDataArr2){
	int timeLength=0;
	int num=0;
	
	FFTObj devFFTObj=NULL;
	int devFFTLength=0;

	float *devDataArr=NULL; // cqt mag&fft mag

	float *devRealArr1=NULL; // fft
	float *devImageArr1=NULL;

	float *devRealArr2=NULL; // ifft
	float *devImageArr2=NULL;

	timeLength=cqtObj->timeLength;
	num=cqtObj->num;

	// 1. 处理缓存
	_cqtObj_dealDeconv(cqtObj);
	devFFTObj=cqtObj->devFFTObj;
	devFFTLength=cqtObj->devFFTLength;

	devDataArr=cqtObj->devDataArr;

	devRealArr1=cqtObj->devRealArr1;
	devImageArr1=cqtObj->devImageArr1;

	devRealArr2=cqtObj->devRealArr2;
	devImageArr2=cqtObj->devImageArr2;

	for(int i=0;i<timeLength;i++){
		// 2. fft&power
		memset(devDataArr, 0, sizeof(float )*devFFTLength);
		memcpy(devDataArr, mDataArr1+i*num, sizeof(float )*num);

		fftObj_fft(devFFTObj, devDataArr, NULL, devRealArr1, devImageArr1);
		__vcabs(devRealArr1, devImageArr1, devFFTLength, devDataArr);

		// 3. timbre --> real(ifft)
		fftObj_ifft(devFFTObj, devDataArr, NULL, devRealArr2, devImageArr2);
		// memcpy(mDataArr2+i*num, devRealArr2, sizeof(float )*num);

		// 4. cqhc
		for(int j=0;j<hcNum;j++){
			int _index=0;

			_index=roundf(cqtObj->binPerOctave*log2f(j+1));
			mDataArr2[i*hcNum+j]=devRealArr2[_index];
		}
	}
}

/***
	mDataArr1 mag/power
	mDataArr2 timbre(formant)
	mDataArr3 pitch
****/
void cqtObj_deconv(CQTObj cqtObj,float *mDataArr1,float *mDataArr2,float *mDataArr3){
	int timeLength=0;
	int num=0;
	
	FFTObj devFFTObj=NULL;
	int devFFTLength=0;

	float *devDataArr=NULL; // cqt mag&fft mag

	float *devRealArr1=NULL; // fft
	float *devImageArr1=NULL;

	float *devRealArr2=NULL; // ifft
	float *devImageArr2=NULL;

	timeLength=cqtObj->timeLength;
	num=cqtObj->num;

	// 1. 处理缓存
	_cqtObj_dealDeconv(cqtObj);
	devFFTObj=cqtObj->devFFTObj;
	devFFTLength=cqtObj->devFFTLength;

	devDataArr=cqtObj->devDataArr;

	devRealArr1=cqtObj->devRealArr1;
	devImageArr1=cqtObj->devImageArr1;

	devRealArr2=cqtObj->devRealArr2;
	devImageArr2=cqtObj->devImageArr2;

	for(int i=0;i<timeLength;i++){
		// 2. fft&mag
		memset(devDataArr, 0, sizeof(float )*devFFTLength);
		memcpy(devDataArr, mDataArr1+i*num, sizeof(float )*num);

		fftObj_fft(devFFTObj, devDataArr, NULL, devRealArr1, devImageArr1);
		__vcabs(devRealArr1, devImageArr1, devFFTLength, devDataArr);

		// 3. timbre --> real(ifft)
		fftObj_ifft(devFFTObj, devDataArr, NULL, devRealArr2, devImageArr2);
		memcpy(mDataArr2+i*num, devRealArr2, sizeof(float )*num);

		// 4. pitch --> real(ifft)
		for(int j=0;j<devFFTLength;j++){
			float _value=0;

			_value=devDataArr[j];
			if(_value<1e-16){
				_value=1e-16;
			}

			// devRealArr1[j]/=(devDataArr[j]+1e-16);
			// devImageArr1[j]/=(devDataArr[j]+1e-16);
			devRealArr1[j]/=_value;
			devImageArr1[j]/=_value;
		}

		fftObj_ifft(devFFTObj, devRealArr1, devImageArr1, devRealArr2, devImageArr2);
		memcpy(mDataArr3+i*num, devRealArr2, sizeof(float )*num);
	}
}

static void _cqtObj_dealDeconv(CQTObj cqtObj){
	FFTObj devFFTObj=NULL;
	int devFFTLength=0;

	float *devDataArr=NULL;

	float *devRealArr1=NULL;
	float *devImageArr1=NULL;

	float *devRealArr2=NULL;
	float *devImageArr2=NULL;

	int num=0;
	int radix2Exp=0;

	num=cqtObj->num;
	devFFTLength=util_ceilPowerTwo(2*num);
	if(devFFTLength!=cqtObj->devFFTLength){
		devFFTObj=cqtObj->devFFTObj;

		devDataArr=cqtObj->devDataArr;

		devRealArr1=cqtObj->devRealArr1;
		devImageArr1=cqtObj->devImageArr1;

		devRealArr2=cqtObj->devRealArr2;
		devImageArr2=cqtObj->devImageArr2;

		fftObj_free(devFFTObj);

		free(devDataArr);

		free(devRealArr1);
		free(devImageArr1);

		free(devRealArr2);
		free(devImageArr2);

		radix2Exp=util_powerTwoBit(devFFTLength);
		fftObj_new(&devFFTObj, radix2Exp);

		devDataArr=__vnew(devFFTLength, NULL);

		devRealArr1=__vnew(devFFTLength, NULL);
		devImageArr1=__vnew(devFFTLength, NULL);

		devRealArr2=__vnew(devFFTLength, NULL);
		devImageArr2=__vnew(devFFTLength, NULL);

		cqtObj->devFFTObj=devFFTObj;
		cqtObj->devFFTLength=devFFTLength;

		cqtObj->devDataArr=devDataArr;

		cqtObj->devRealArr1=devRealArr1;
		cqtObj->devImageArr1=devImageArr1;

		cqtObj->devRealArr2=devRealArr2;
		cqtObj->devImageArr2=devImageArr2;
	}
}

static void _cqtObj_cqt(CQTObj cqtObj,float *dataArr,int dataLength,float *mRealArr3,float *mImageArr3){
	STFTObj stftObj=NULL;
	ResampleObj resampleObj=NULL;

	int fftLength=0; // fftLength,timeLength,num
	int timeLength=0;
	int stftTimeLength=0;

	int slideLength=0;

	int num=0;
	int binPerOctave=0;
	int octaveNum=0;

	float *mRealFilterBankArr=NULL; // num*(fftLength/2+1)
	float *mImageFilterBankArr=NULL;

	float *freBandArr=NULL; // num
	float *lenArr=NULL;

	float *sLenArr=NULL;
	float *dLenArr=NULL;
	
	float *mRealArr1=NULL;
	float *mImageArr1=NULL;

	float *mRealArr2=NULL;
	float *mImageArr2=NULL;

	int normType=0;

	float *preDataArr=NULL;
	float *curDataArr=NULL;

	int preDataLength=0;
	int curDataLength=0;

	int index1=0;
	int vFlag=0;

	int isContinue=0;

	vFlag=cqtObj->vFlag;
	isContinue=cqtObj->isContinue;

	stftObj=cqtObj->stftObj;
	resampleObj=cqtObj->resampleObj;

	fftLength=cqtObj->fftLength;
	
	num=cqtObj->num;
	binPerOctave=cqtObj->binPerOctave;
	octaveNum=cqtObj->octaveNum;

	mRealFilterBankArr=cqtObj->mRealFilterBankArr; // num*(fftLength/2+1)
	mImageFilterBankArr=cqtObj->mImageFilterBankArr;

	freBandArr=cqtObj->freBandArr;
	lenArr=cqtObj->lenArr;

	sLenArr=cqtObj->sLenArr;
	dLenArr=cqtObj->dLenArr;

	mRealArr1=cqtObj->mRealArr1; // timeLength*fftLength
	mImageArr1=cqtObj->mImageArr1;

	mRealArr2=cqtObj->mRealArr2; // timeLength*binPerOctave
	mImageArr2=cqtObj->mImageArr2;

	normType=cqtObj->normType;

	preDataArr=__vnew(dataLength, NULL);
	curDataArr=__vnew(dataLength, NULL);
	preDataLength=dataLength;

	slideLength=cqtObj->slideLength;

	stftTimeLength=dataLength/slideLength+1;
	if(!isContinue){
		timeLength=stftTimeLength;
	}
	else{
		timeLength=(dataLength-fftLength)/slideLength+1;
	}

	if(cqtObj->stftTimeLength<stftTimeLength||
		cqtObj->stftTimeLength>stftTimeLength*2){ // 更新缓存
		free(mRealArr1);
		free(mImageArr1);

		free(mRealArr2);
		free(mImageArr2);

		/***
			down stft timeLength可能比top stft timeLength多1
			cqt timeLength以top stft timeLength为准
			内存尺度 timeLength*fftLength => (timeLength+1)*fftLength
		****/
		mRealArr1=__vnew((stftTimeLength+1)*fftLength, NULL);
		mImageArr1=__vnew((stftTimeLength+1)*fftLength, NULL);

		mRealArr2=__vnew((stftTimeLength+1)*binPerOctave, NULL);
		mImageArr2=__vnew((stftTimeLength+1)*binPerOctave, NULL);
	}

	// 1. top octave
	stftObj_setSlideLength(stftObj, slideLength);
	stftObj_stft(stftObj,dataArr,dataLength,mRealArr1,mImageArr1);

	// t*fftLength => t*(fftLength/2+1)
	for(int i=1;i<timeLength;i++){
		for(int j=0;j<(fftLength/2+1);j++){
			mRealArr1[i*(fftLength/2+1)+j]=mRealArr1[i*fftLength+j];
			mImageArr1[i*(fftLength/2+1)+j]=mImageArr1[i*fftLength+j];
		}
	}

	// S.dot(filterBank) t*(fftLength/2+1) @ num*(fftLength/2+1)
	if(vFlag){ // vqt num
		index1=(octaveNum-1)*binPerOctave*(fftLength/2+1);
	}
	__mcdot1(mRealArr1, mImageArr1, 
			mRealFilterBankArr+index1, mImageFilterBankArr+index1,
		 	timeLength, fftLength/2+1,
		 	binPerOctave, fftLength/2+1, 
		 	mRealArr2, mImageArr2);

	for(int i=0;i<timeLength;i++){
		for(int j=(octaveNum-1)*binPerOctave,k=0;j<octaveNum*binPerOctave;j++,k++){
			float _value1=0;
			float _value2=0;

			_value1=mRealArr2[i*binPerOctave+k];
			_value2=mImageArr2[i*binPerOctave+k];

			// norm => similar fft 'orthi' sqrt(len)
			if(cqtObj->isScale){
				_value1/=sLenArr[j];
				_value2/=sLenArr[j];
			}

			mRealArr3[i*octaveNum*binPerOctave+j]=_value1;
			mImageArr3[i*octaveNum*binPerOctave+j]=_value2;
		}
	}

	// 2. down
	memcpy(preDataArr, dataArr, sizeof(float )*dataLength);
	for(int i=octaveNum-2;i>=0;i--){
		// resample
		curDataLength=resampleObj_resample(resampleObj, preDataArr, preDataLength, curDataArr);
		slideLength/=2;

		// stft t*fftLength
		stftObj_setSlideLength(stftObj, slideLength);
		stftObj_stft(stftObj,curDataArr,curDataLength,mRealArr1,mImageArr1);

		// t*fftLength => t*(fftLength/2+1)
		for(int j=1;j<timeLength;j++){
			for(int k=0;k<(fftLength/2+1);k++){
				mRealArr1[j*(fftLength/2+1)+k]=mRealArr1[j*fftLength+k];
				mImageArr1[j*(fftLength/2+1)+k]=mImageArr1[j*fftLength+k];
			}
		}

		// S.dot(filterBank) t*(fftLength/2+1) @ num*(fftLength/2+1)
		if(vFlag){ // vqt num
			index1=i*binPerOctave*(fftLength/2+1);
		}
		__mcdot1(mRealArr1, mImageArr1, 
				mRealFilterBankArr+index1, mImageFilterBankArr+index1,
			 	timeLength, fftLength/2+1,
			 	binPerOctave, fftLength/2+1, 
			 	mRealArr2, mImageArr2);

		for(int n=0;n<timeLength;n++){
			for(int j=i*binPerOctave,k=0;j<(i+1)*binPerOctave;j++,k++){
				float _value1=0;
				float _value2=0;

				_value1=mRealArr2[n*binPerOctave+k];
				_value2=mImageArr2[n*binPerOctave+k];

				// norm => downsample 默认非scale sqrt(2^i)
				_value1*=dLenArr[octaveNum-i-1];
				_value2*=dLenArr[octaveNum-i-1];

				// norm => similar fft 'orthi' sqrt(len)
				if(cqtObj->isScale){
					_value1/=sLenArr[j];
					_value2/=sLenArr[j];
				}

				mRealArr3[n*octaveNum*binPerOctave+j]=_value1;
				mImageArr3[n*octaveNum*binPerOctave+j]=_value2;
			}
		}

		// deal data
		memset(preDataArr, 0, sizeof(float )*preDataLength);
		memcpy(preDataArr, curDataArr, sizeof(float )*curDataLength);
		memset(curDataArr, 0, sizeof(float )*curDataLength);
		preDataLength=curDataLength;
	}

	free(preDataArr);
	free(curDataArr);

	cqtObj->timeLength=timeLength;
	cqtObj->stftTimeLength=stftTimeLength;

	cqtObj->mRealArr1=mRealArr1;
	cqtObj->mImageArr1=mImageArr1;

	cqtObj->mRealArr2=mRealArr2;
	cqtObj->mImageArr2=mImageArr2;
}

void cqtObj_free(CQTObj cqtObj){
	STFTObj stftObj=NULL; // octaveNum
	ResampleObj resampleObj=NULL;

	float *mRealFilterBankArr=NULL; // num*(fftLength/2+1)
	float *mImageFilterBankArr=NULL;

	float *freBandArr=NULL; // num
	float *lenArr=NULL; // binPerOctave

	float *sLenArr=NULL;
	float *dLenArr=NULL;

	float *mRealArr1=NULL; // stft r,i(timeLength*fftLength) -> power r(timeLength*(fftLength/2+1)) 
	float *mImageArr1=NULL;

	float *mRealArr2=NULL; 
	float *mImageArr2=NULL;

	float *tailDataArr=NULL;
	float *validDataArr=NULL;

	float *mChromaFilterBank=NULL;
	float *mSArr=NULL;

	FFTObj fftObj=NULL; // dct加速
	DCTObj dctObj=NULL;

	FFTObj devFFTObj=NULL;

	float *devDataArr=NULL; // spectral mag&fft mag

	float *devRealArr1=NULL; // fft
	float *devImageArr1=NULL;

	float *devRealArr2=NULL; // ifft
	float *devImageArr2=NULL;

	if(!cqtObj){
		return;
	}

	stftObj=cqtObj->stftObj;
	resampleObj=cqtObj->resampleObj;

	mRealFilterBankArr=cqtObj->mRealFilterBankArr;
	mImageFilterBankArr=cqtObj->mImageFilterBankArr;

	freBandArr=cqtObj->freBandArr;
	lenArr=cqtObj->lenArr;

	sLenArr=cqtObj->sLenArr;
	dLenArr=cqtObj->dLenArr;

	mRealArr1=cqtObj->mRealArr1;
	mImageArr1=cqtObj->mImageArr1;

	mRealArr2=cqtObj->mRealArr2;
	mImageArr2=cqtObj->mImageArr2;

	tailDataArr=cqtObj->tailDataArr;
	validDataArr=cqtObj->validDataArr;

	mChromaFilterBank=cqtObj->mChromaFilterBank;
	mSArr=cqtObj->mSArr;

	fftObj=cqtObj->fftObj;
	dctObj=cqtObj->dctObj;

	devFFTObj=cqtObj->devFFTObj;
	devDataArr=cqtObj->devDataArr;

	devRealArr1=cqtObj->devRealArr1;
	devImageArr1=cqtObj->devImageArr1;

	devRealArr2=cqtObj->devRealArr2;
	devImageArr2=cqtObj->devImageArr2;

	stftObj_free(stftObj);
	resampleObj_free(resampleObj);

	free(mRealFilterBankArr);
	free(mImageFilterBankArr);

	free(freBandArr);
	free(lenArr);

	free(sLenArr);
	free(dLenArr);

	free(mRealArr1);
	free(mImageArr1);

	free(mRealArr2);
	free(mImageArr2);

	free(tailDataArr);
	free(validDataArr);

	free(mChromaFilterBank);
	free(mSArr);

	fftObj_free(fftObj);
	dctObj_free(dctObj);

	fftObj_free(devFFTObj);
	free(devDataArr);

	free(devRealArr1);
	free(devImageArr1);

	free(devRealArr2);
	free(devImageArr2);

	free(cqtObj);
}

// filterBank相关数据
static void _cqtObj_dealFilterBank(CQTObj cqtObj,int num,float minFre,int samplate,
								int binPerOctave,SpectralFilterBankNormalType normType,WindowType winType,
								float factor,float beta,float thresh,int vFlag){
	int fftLength=0; 
	
	float *freBandArr=NULL;
	float *lenArr=NULL;

	float *sLenArr=NULL;
	float *dLenArr=NULL;

	float *mRealFilterBankArr=NULL; // num*(fftLength/2+1)
	float *mImageFilterBankArr=NULL;

	int octaveNum=0;
	int index=0;

	octaveNum=num/binPerOctave;
	index=(octaveNum-1)*binPerOctave;

	lenArr=__vnew(binPerOctave, NULL);
	sLenArr=__vnew(num, NULL);
	dLenArr=__vnew(octaveNum, NULL);

	freBandArr=cqt_calFreArr(minFre,num,binPerOctave);
	fftLength=cqt_calFFTLength(freBandArr[index],samplate,binPerOctave,&factor,&beta);

	// vFlag 1 => vqt 设num 
	mRealFilterBankArr=__vnew((vFlag?num:binPerOctave)*(fftLength/2+1), NULL);
	mImageFilterBankArr=__vnew((vFlag?num:binPerOctave)*(fftLength/2+1), NULL);
	
	cqt_calLengthArr(binPerOctave,freBandArr+index,samplate,
					binPerOctave,
					&factor,&beta,
					lenArr);

	cqt_calLengthArr(num,freBandArr,samplate,
					binPerOctave,
					&factor,&beta,
					sLenArr);

	for(int i=0;i<num;i++){
		sLenArr[i]=sqrtf(sLenArr[i]);
	}

	dLenArr[0]=1;
	for(int i=1;i<octaveNum;i++){
		dLenArr[i]=sqrtf(1<<i);
	}

	// ??? cqt_filterBank替换
	// cqt_downFilterBank(num,freBandArr,samplate,
	// 				binPerOctave,normType,&winType,
	// 				&factor,&beta,&thresh,
	// 				lenArr,fftLength,
	// 				mRealFilterBankArr,mImageFilterBankArr);

	// vFlag 1 => vqt 设num 不能共用top octave freKernal
	index=(vFlag?0:index);
	cqt_downFilterBank(vFlag?num:binPerOctave,freBandArr+index,samplate,
					binPerOctave,normType,&winType,
					&factor,&beta,&thresh,
					lenArr,fftLength,
					mRealFilterBankArr,mImageFilterBankArr);

	cqtObj->fftLength=fftLength;

	cqtObj->num=num;
	cqtObj->octaveNum=octaveNum;
	cqtObj->binPerOctave=binPerOctave;

	cqtObj->mRealFilterBankArr=mRealFilterBankArr;
	cqtObj->mImageFilterBankArr=mImageFilterBankArr;

	cqtObj->freBandArr=freBandArr;
	cqtObj->lenArr=lenArr;

	cqtObj->sLenArr=sLenArr;
	cqtObj->dLenArr=dLenArr;

	cqtObj->samplate=samplate;

	cqtObj->windowType=winType;
	cqtObj->normType=normType;
}

static void _cqtObj_dealResample(CQTObj cqtObj){
	ResampleObj resampleObj=NULL;

	ResampleQualityType type=ResampleQuality_Fast;
	int isScale=1;

	resampleObj_new(&resampleObj,&type,&isScale,NULL);
	resampleObj_setSamplate(resampleObj,2,1); // 降采样

	cqtObj->resampleObj=resampleObj;
}

// stft数据相关 winType rect
static void _cqtObj_dealStftArr(CQTObj cqtObj,int fftLength,int slideLength){
	STFTObj *stftObjArr=NULL;

	int octaveNum=0;
	int radix2Exp=0;

	octaveNum=cqtObj->octaveNum;
	radix2Exp=util_powerTwoBit(fftLength);

	cqtObj->radix2Exp=radix2Exp;
	cqtObj->slideLength=slideLength;
	
	stftObjArr=(STFTObj *)calloc(octaveNum, sizeof(STFTObj ));
	stftObj_new(stftObjArr,radix2Exp,NULL,&slideLength,NULL);

	for(int i=1;i<octaveNum;i++){
		slideLength/=2;
		stftObj_new(stftObjArr+i,radix2Exp,NULL,&slideLength,NULL);
	}	

	// cqtObj->stftObjArr=stftObjArr;
}

static void _cqtObj_dealStft(CQTObj cqtObj,int fftLength,int slideLength,int isContinue){
	STFTObj stftObj=NULL;

	PaddingPositionType pType=PaddingPosition_Right;
	int radix2Exp=0;

	radix2Exp=util_powerTwoBit(fftLength);

	cqtObj->radix2Exp=radix2Exp;
	cqtObj->slideLength=slideLength;
	
	stftObj_new(&stftObj,radix2Exp,NULL,&slideLength,NULL);
	stftObj_enablePadding(stftObj,1); // 默认 center 0 padding

	if(isContinue){ // right 0 padding
		stftObj_setPadding(stftObj, &pType, NULL, NULL, NULL);
	}

	cqtObj->stftObj=stftObj;
}

static int _calRadix2(int length){
	int r=0;
	int m=0;

	while(1){
		m=length%2;
		if(!m){ //
			r++;
			length=length/2;
			if(length==1){
				break;
			}
		}
		else{
			r=0;
			break;
		}
	}

	return r;
}



