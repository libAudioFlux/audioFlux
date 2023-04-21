// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/flux_correct.h"
#include "dsp/fft_algorithm.h"
#include "dsp/dct_algorithm.h"

#include "filterbank/auditory_weight.h"
#include "filterbank/auditory_filterBank.h"
#include "filterbank/chroma_filterBank.h"

#include "flux_spectral.h"
#include "stft_algorithm.h"
#include "spectrogram_algorithm.h"

typedef struct {
	float maxMin;
	float minMax;

	float ratio;

	int isRecover;
	int isKWeight;

} SpectralDeepConfig;

struct OpaqueSpectrogram{
	STFTObj stftObj;

	int fftLength; // fftLength,timeLength,num
	int timeLength;

	int num;
	float *mFilterBankArr; // num*(fftLength/2+1)/baseNum*(fftLength/2+1)
	float *mChromaFilterBankArr; // num*baseNum

	float *freBandArr;
	int *binBandArr;

	// stft相关
	float *mRealArr; // stft r,i(timeLength*fftLength) -> power r(timeLength*(fftLength/2+1)) 
	float *mImageArr;
	float *mSArr; // timeLength*(fftLength/2+1)

	float *energyArr; // timeLength
	float *vArr1; // fftLength

	int *indexArr;
	int indexLength;

	int start; // 默认0~num-1 可选区间计算spectral 0~baseNum-1???
	int end;

	float *sumArr; // timeLength 内存和值都缓存
	float *cArr1; // μ1
	float *cArr2;
	float *entropyArr;

	float *meanFreArr; // 仅内存缓存
	float *meanValueArr; 
	
	// 重计算标志 1. 内存刷新 2. start-end变换 3. data发生变化呢
	int isSum; 
	int isC1;
	int isC2;
	int isEntropy; // isEntropy&&isEnNorm 
	int isEnNorm;

	int isMean;

	int isDebug;

	// dct相关
	FFTObj fftObj; // dct加速
	DCTObj dctObj; // 直接矩阵

	// params相关
	int samplate; // filterBank 
	float lowFre;
	float highFre;

	int lowIndex; // 针对linear
	int highIndex;

	int binPerOctave; // 针对Log/LogChroma
	/***
		1. Chroma baseNum-->freBandArr
		1. LogChroma baseFre-->chrome filterBank/freBandArr
		2. Deep baseNum&baseFre-->indexArr/correctFreArr
		3. DeepChroma +chroma filterBank
	****/
	int baseNum; 
	float baseFre;

	int radix2Exp; // stft
	WindowType windowType;
	int slideLength;
	int isContinue;

	SpectralDataType dataType; // data&&filterBank
	SpectralFilterBankScaleType filterScaleType;
	SpectralFilterBankStyleType filterStyleType;
	SpectralFilterBankNormalType filterNormalType;
	ChromaDataNormalType dataNormType; // 针对stft-chroma

	float normValue;  // 默认1 >0 mag=>result power=>S

	// deconv相关
	FFTObj devFFTObj;
	int devFFTLength;

	float *devDataArr; // spectral mag&fft mag

	float *devRealArr1; // fft
	float *devImageArr1;

	float *devRealArr2; // ifft
	float *devImageArr2;

	// deep spectrogram相关 PCP/HPCP 基于amp的过滤模型
	float *ampArr; // fftLength/2+1
	float *weightArr;

	int *salienceIndexArr; //  针对deep/deepChroma
	int salienceLength; // max==endIndex-startIndex+1

	int startIndex; // stft index
	int endIndex;

	int midiStart; // lowFre~highFre start-end+1==baseNum
	int midiEnd;

	int deepOrder; // 1

	float maxMin; // 13.0
	float minMax; // 2.0

	float ratio; // 10.0

	int isRecover; // 0 是否恢复amp
	int isKWeight; // 0 是否K计权

	float *mCorrectFreArr; // 校正频率 timeLength*baseNum
	float *mToneFreArr;
	float *mToneFreArr1; // 临近参照
};

void spectrogramObj_setDeepConfig(SpectrogramObj spectrogramObj,SpectralDeepConfig *config);

static void __spectrogramObj_init(SpectrogramObj spectrogramObj);
// 计算高频sum/c1/c2
static void __spectrogramObj_calSum(SpectrogramObj spectrogramObj,float *mDataArr);
static void __spectrogramObj_calC1(SpectrogramObj spectrogramObj,float *mDataArr);
static void __spectrogramObj_calC2(SpectrogramObj spectrogramObj,float *mDataArr);
static void __spectrogramObj_calEntropy(SpectrogramObj spectrogramObj,float *mDataArr,int isNorm);
static void __spectrogramObj_calMean(SpectrogramObj spectrogramObj,float *mDataArr);

static void __spectrogramObj_spectrogram(SpectrogramObj spectrogramObj,float *dataArr,int dataLength,
									float *mRealArr,float *mImageArr,int nLength,int mLength,
									float *mSpectArr,float *mPhaseArr);

// mfcc/gtcc ???
static void __spectrogramObj_xxcc(SpectrogramObj spectrogramObj,float *mDataArr1,int mLength,CepstralRectifyType *rectifyType,float *mDataArr2);
static void __spectrogramObj_dealDeconv(SpectrogramObj spectrogramObj);
static void __spectrogramObj_deepFilter(SpectrogramObj spectrogramOb,float *mDataArr1,float *mDataArr2,int isDeep);

static int __spectrogramObj_calSalience(SpectrogramObj spectrogramObj,float *ampDataArr);
static void __spectrogramObj_calLinearBandArr(SpectrogramObj spectrogramObj,float *freBandArr,int *binBandArr);
static void __spectrogramObj_calLogBandArr(SpectrogramObj spectrogramObj,float *freBandArr,int *binBandArr);
static void __spectrogramObj_calDeepBandArr(SpectrogramObj spectrogramObj,float *freBandArr,int *binBandArr);

static void _calKWeight(float *ampArr1,int start,int end,int fftLength,float *weightArr,float *ampArr2);

static int _calBaseNum(float lowFre,float highFre,float binPerOctave);
static float _calBaseFre(float lowFre,float binPerOctave);
static void _calTone(float value1,float *value2,float *value3);
static int _calRadix2(int length);

int spectrogramObj_newLinear(SpectrogramObj *spectrogramObj,int samplate,int radix2Exp,int *isContinue){
	int status=0;
	
	int _samplate=0;
	int _radix2Exp=0;
	SpectralFilterBankScaleType filterScaleType;

	_samplate=samplate;
	_radix2Exp=radix2Exp;
	filterScaleType=SpectralFilterBankScale_Linear;
	
	status=spectrogramObj_new(spectrogramObj,2,
							&_samplate,NULL,NULL,NULL,
							&_radix2Exp,NULL,NULL,isContinue,
							NULL,&filterScaleType,NULL,NULL);

	return status;
}

int spectrogramObj_newMel(SpectrogramObj *spectrogramObj,int num,int samplate,int radix2Exp,int *isContinue){
	int status=0;

	int _samplate=0;
	int _radix2Exp=0;
	SpectralFilterBankScaleType filterScaleType;

	_samplate=samplate;
	_radix2Exp=radix2Exp;
	filterScaleType=SpectralFilterBankScale_Mel;
	
	status=spectrogramObj_new(spectrogramObj,num,
							&_samplate,NULL,NULL,NULL,
							&_radix2Exp,NULL,NULL,isContinue,
							NULL,&filterScaleType,NULL,NULL);

	return status;
}

int spectrogramObj_newBark(SpectrogramObj *spectrogramObj,int num,int samplate,int radix2Exp,int *isContinue){
	int status=0;
	
	int _samplate=0;
	int _radix2Exp=0;
	SpectralFilterBankScaleType filterScaleType;

	_samplate=samplate;
	_radix2Exp=radix2Exp;
	filterScaleType=SpectralFilterBankScale_Bark;
	
	status=spectrogramObj_new(spectrogramObj,num,
							&_samplate,NULL,NULL,NULL,
							&_radix2Exp,NULL,NULL,isContinue,
							NULL,&filterScaleType,NULL,NULL);


	return status;
}

int spectrogramObj_newErb(SpectrogramObj *spectrogramObj,int num,int samplate,int radix2Exp,int *isContinue){
	int status=0;
	
	int _samplate=0;
	int _radix2Exp=0;
	SpectralFilterBankScaleType filterScaleType;

	_samplate=samplate;
	_radix2Exp=radix2Exp;
	filterScaleType=SpectralFilterBankScale_Erb;
	
	status=spectrogramObj_new(spectrogramObj,num,
							&_samplate,NULL,NULL,NULL,
							&_radix2Exp,NULL,NULL,isContinue,
							NULL,&filterScaleType,NULL,NULL);


	return status;
}

int spectrogramObj_newChroma(SpectrogramObj *spectrogramObj,int samplate,int radix2Exp,int *isContinue){
	int status=0;
	
	int num=12;
	int _samplate=0;
	int _radix2Exp=0;
	SpectralFilterBankScaleType filterScaleType;

	_samplate=samplate;
	_radix2Exp=radix2Exp;
	filterScaleType=SpectralFilterBankScale_Chroma;
	
	status=spectrogramObj_new(spectrogramObj,num,
							&_samplate,NULL,NULL,NULL,
							&_radix2Exp,NULL,NULL,isContinue,
							NULL,&filterScaleType,NULL,NULL);


	return status;
}

int spectrogramObj_newDeep(SpectrogramObj *spectrogramObj,int num,int samplate,int radix2Exp,int *isContinue){
	int status=0;
	
	int _samplate=0;
	int _radix2Exp=0;
	SpectralFilterBankScaleType filterScaleType;

	_samplate=samplate;
	_radix2Exp=radix2Exp;
	filterScaleType=SpectralFilterBankScale_Deep;
	
	status=spectrogramObj_new(spectrogramObj,num,
							&_samplate,NULL,NULL,NULL,
							&_radix2Exp,NULL,NULL,isContinue,
							NULL,&filterScaleType,NULL,NULL);


	return status;
}

int spectrogramObj_newDeepChroma(SpectrogramObj *spectrogramObj,int samplate,int radix2Exp,int *isContinue){
	int status=0;
	
	int num=12;
	int _samplate=0;
	int _radix2Exp=0;
	SpectralFilterBankScaleType filterScaleType;

	_samplate=samplate;
	_radix2Exp=radix2Exp;
	filterScaleType=SpectralFilterBankScale_DeepChroma;
	
	status=spectrogramObj_new(spectrogramObj,num,
							&_samplate,NULL,NULL,NULL,
							&_radix2Exp,NULL,NULL,isContinue,
							NULL,&filterScaleType,NULL,NULL);


	return status;
}

int spectrogramObj_new(SpectrogramObj *spectrogramObj,int num,
					int *samplate,float *lowFre,float *highFre,int *binPerOctave,
					int *radix2Exp,WindowType *windowType,int *slideLength,int *isContinue,
					SpectralDataType *dataType,
					SpectralFilterBankScaleType *filterScaleType,
					SpectralFilterBankStyleType *filterStyleType,
					SpectralFilterBankNormalType *filterNormalType){
	int status=0;
	SpectrogramObj spec=NULL;

	int fftLength=0;

	int _samplate=32000; // filterBank 
	float _lowFre=0;
	float _highFre=0;

	int lowIndex=0; // 针对linear
	int highIndex=0;

	int midiStart=0; // 针对deep/deepChroma
	int midiEnd=0;

	int _binPerOctave=12; // 针对log

	int _radix2Exp=12; // stft
	WindowType _windowType=Window_Hann;
	int _slideLength=0;
	int _isContinue=0;

	SpectralDataType _dataType=SpectralData_Power; // data&&filterBank
	SpectralFilterBankScaleType _filterScaleType=SpectralFilterBankScale_Linear;
	SpectralFilterBankStyleType _filterStyleType=SpectralFilterBankStyle_Slaney;
	SpectralFilterBankNormalType _filterNormalType=SpectralFilterBankNormal_None;

	int baseNum=0;
	float baseFre=0; // 针对log-chroma/deep-chroma 

	if(radix2Exp){
		_radix2Exp=*radix2Exp;
		if(_radix2Exp<1||_radix2Exp>30){
			status=-100;
			printf("radix2Exp is error!\n");
			return status;
		}
	}

	fftLength=1<<_radix2Exp;
	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
		}
	}

	if(dataType){
		_dataType=*dataType;
	}
	if(filterScaleType){
		_filterScaleType=*filterScaleType;
	}
	if(filterStyleType){
		_filterStyleType=*filterStyleType;
	}
	if(filterNormalType){
		_filterNormalType=*filterNormalType;
	}

	_highFre=_samplate/2.0;

	if(lowFre){
		if(*lowFre>=0&&*lowFre<_samplate/2.0){
			_lowFre=*lowFre;
		}
	}

	if(_lowFre==0){
		if(_filterScaleType==SpectralFilterBankScale_Octave||
			_filterScaleType==SpectralFilterBankScale_Log||
			_filterScaleType==SpectralFilterBankScale_LogChroma||
			_filterScaleType==SpectralFilterBankScale_Deep||
			_filterScaleType==SpectralFilterBankScale_DeepChroma){ // Log/LogChroma/Deep/DeepChroma

			_lowFre=powf(2, -45/12.0)*440;
			_highFre=powf(2, 38/12.0)*440;
		}
	}

	if(highFre){
		if(*highFre>0&&*highFre<=_samplate/2.0){
			_highFre=*highFre;
		}
	}

	if(_highFre<_lowFre){
		_lowFre=0;
		_highFre=_samplate/2.0;
		if(_filterScaleType==SpectralFilterBankScale_Octave||
			_filterScaleType==SpectralFilterBankScale_Log||
			_filterScaleType==SpectralFilterBankScale_LogChroma||
			_filterScaleType==SpectralFilterBankScale_Deep||
			_filterScaleType==SpectralFilterBankScale_DeepChroma){ // Log/LogChroma/Deep/DeepChroma

			_lowFre=powf(2, -45/12.0)*440;
			_highFre=powf(2, 38/12.0)*440;
		}
	}

	if(binPerOctave){
		if(*binPerOctave>0){
			_binPerOctave=*binPerOctave;
		}
	}

	if(_binPerOctave%12!=0){
		_binPerOctave=12;
	}

	if(_filterScaleType==SpectralFilterBankScale_Linear||
		_filterScaleType==SpectralFilterBankScale_Chroma){
		float det=0;

		det=_samplate/(float )fftLength;
		lowIndex=roundf(_lowFre/det);
		highIndex=roundf(_highFre/det);
	}

	if(_filterScaleType==SpectralFilterBankScale_Deep||
		_filterScaleType==SpectralFilterBankScale_DeepChroma){ // deep/deepChroma hamm
		_windowType=Window_Hamm;
	}

	if(windowType){
		_windowType=*windowType;
	}

	if(_filterScaleType==SpectralFilterBankScale_Deep||
		_filterScaleType==SpectralFilterBankScale_DeepChroma){
		if(_windowType>Window_Hamm){
			_windowType=Window_Hamm;
		}
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

	if(_filterScaleType==SpectralFilterBankScale_Linear){ // num 和low/high有关 ???
		num=highIndex-lowIndex+1;
	}
	else if(_filterScaleType==SpectralFilterBankScale_Octave){ // Log
		auditory_reviseLogFre(num, _lowFre,_highFre, _binPerOctave,1, &_lowFre, &_highFre);

		if(_highFre>_samplate/2.0){
			printf("scale log: lowFre and num is large, overflow error!\n");
			return -1;
		}

		baseNum=num;
		baseFre=_lowFre;
	}
	else if(_filterScaleType==SpectralFilterBankScale_Deep){ // deep
		baseNum=num;
		baseFre=_calBaseFre(_lowFre, 12);

		midiStart=auditory_freToMidi(baseFre);
		midiEnd=baseNum+midiStart-1;
	}
	else if(_filterScaleType==SpectralFilterBankScale_Chroma){ // stft-chroma处理 不需要binPerOctave倍数关系
		if(num<12||num%12!=0){  // 12/24/36/48...
			num=12;
		}

		baseNum=highIndex-lowIndex+1;
	}
	else if(_filterScaleType==SpectralFilterBankScale_LogChroma){ // log-chroma处理
		if(num<=0){
			num=12;
		}
		else if(num>_binPerOctave||_binPerOctave%num!=0){
			num=12;
		}

		baseNum=_calBaseNum(_lowFre, _highFre, _binPerOctave);
		baseFre=_calBaseFre(_lowFre, _binPerOctave);
	}
	else if(_filterScaleType==SpectralFilterBankScale_DeepChroma){ // DeepChroma
		if(num<12||num%12!=0){
			num=12;
		}

		baseNum=_calBaseNum(_lowFre, _highFre, 12);
		baseFre=_calBaseFre(_lowFre, 12);

		midiStart=auditory_freToMidi(baseFre);
		midiEnd=baseNum+midiStart-1;
	}

	if(num<2||num>fftLength/2+1){ 
		printf("num is error!\n");
		return -1;
	}

	spec=*spectrogramObj=(SpectrogramObj )calloc(1, sizeof(struct OpaqueSpectrogram ));

	spec->fftLength=fftLength;
	spec->num=num;

	spec->samplate=_samplate;
	spec->lowFre=_lowFre;
	spec->highFre=_highFre;

	spec->lowIndex=lowIndex;
	spec->highIndex=highIndex;

	spec->binPerOctave=_binPerOctave;
	spec->baseNum=baseNum;
	spec->baseFre=baseFre;

	spec->radix2Exp=_radix2Exp;
	spec->windowType=_windowType;
	spec->slideLength=_slideLength;
	spec->isContinue=_isContinue;

	spec->dataType=_dataType;
	spec->filterScaleType=_filterScaleType;
	spec->filterStyleType=_filterStyleType;
	spec->filterNormalType=_filterNormalType;
	spec->dataNormType=ChromaDataNormal_Max;

	spec->normValue=1;

	// deep相关
	spec->deepOrder=1;

	spec->midiStart=midiStart;
	spec->midiEnd=midiEnd;

	spec->maxMin=13.0;
	spec->minMax=2.0;

	spec->ratio=10.0;

	__spectrogramObj_init(spec);

	return status;
}

/***
	1. stft
	2. filterBank/chromaFilterBank
	3. dct
	4. deep
****/
static void __spectrogramObj_init(SpectrogramObj spectrogramObj){
	STFTObj stftObj=NULL;

	int num=0;
	int fftLength=0;

	float *mFilterBankArr=NULL;
	float *mChromaFilterBankArr=NULL;

	float *freBandArr=NULL;
	int *binBandArr=NULL;

	float *vArr1=NULL;

	int start=0;
	int end=0;

	int samplate=0; // filterBank 

	float lowFre=0;
	float highFre=0;

	int lowIndex=0; 
	int highIndex=0; 

	int binPerOctave=0;
	int baseNum=0;
	float baseFre=0;

	int radix2Exp=0; // stft
	WindowType windowType;
	int slideLength=0;
	int isContinue=0;

	FFTObj fftObj=NULL;
	DCTObj dctObj=NULL;
	int r=0;

	float *ampArr=NULL; // fftLength/2+1
	float *weightArr=NULL;

	int *salienceIndexArr=NULL; //  针对deep/deepChroma

	int startIndex=0; 
	int endIndex=0;
	float det=0;

	int *indexArr=NULL;
	int indexLength=0;

	SpectralFilterBankScaleType filterScaleType;
	SpectralFilterBankStyleType filterStyleType;
	SpectralFilterBankNormalType filterNormalType;

	num=spectrogramObj->num;
	fftLength=spectrogramObj->fftLength;

	samplate=spectrogramObj->samplate;
	lowFre=spectrogramObj->lowFre;
	highFre=spectrogramObj->highFre;

	lowIndex=spectrogramObj->lowIndex;
	highIndex=spectrogramObj->highIndex;

	binPerOctave=spectrogramObj->binPerOctave;
	baseNum=spectrogramObj->baseNum;
	baseFre=spectrogramObj->baseFre;

	radix2Exp=spectrogramObj->radix2Exp;
	windowType=spectrogramObj->windowType;
	slideLength=spectrogramObj->slideLength;
	isContinue=spectrogramObj->isContinue;

	filterScaleType=spectrogramObj->filterScaleType;
	filterStyleType=spectrogramObj->filterStyleType;
	filterNormalType=spectrogramObj->filterNormalType;

	vArr1=__vnew(fftLength, NULL);

	// 1. stft
	stftObj_new(&stftObj,radix2Exp,&windowType,&slideLength,&isContinue);

	// 2. filterBank
	if(filterScaleType!=SpectralFilterBankScale_Linear){ // mel/bark/erb... chroma

		if(filterScaleType==SpectralFilterBankScale_Mel||
			filterScaleType==SpectralFilterBankScale_Bark||
			filterScaleType==SpectralFilterBankScale_Erb||
			filterScaleType==SpectralFilterBankScale_Octave||
			filterScaleType==SpectralFilterBankScale_Linspace||
			filterScaleType==SpectralFilterBankScale_Log){ // mel/bark/erb/log linespace/logspace

			mFilterBankArr=__vnew(num*(fftLength/2+1), NULL);
			freBandArr=__vnew(num+2, NULL); 
			binBandArr=__vnewi(num+2, NULL);

			auditory_filterBank(num,fftLength,samplate,0,
								filterScaleType,filterStyleType,filterNormalType,
								lowFre,highFre,binPerOctave,
								mFilterBankArr,
								freBandArr,
								binBandArr);
		}
		else if(filterScaleType==SpectralFilterBankScale_Chroma){ // stft-chroma
			mFilterBankArr=__vnew(num*(fftLength/2+1), NULL);
			freBandArr=__vnew(baseNum+2, NULL); 
			binBandArr=__vnewi(baseNum+2, NULL);

			chroma_stftFilterBank(num,fftLength,samplate,
								NULL,NULL,
								mFilterBankArr);

			__spectrogramObj_calLinearBandArr(spectrogramObj,freBandArr,binBandArr);
		}
		else if(filterScaleType==SpectralFilterBankScale_LogChroma){ // Log-Chroma
			mFilterBankArr=__vnew(baseNum*(fftLength/2+1), NULL);
			freBandArr=__vnew(baseNum+2, NULL); 
			binBandArr=__vnewi(baseNum+2, NULL);

			auditory_filterBank(baseNum,fftLength,samplate,0,
								filterScaleType,filterStyleType,filterNormalType,
								lowFre,highFre,binPerOctave,
								mFilterBankArr,
								freBandArr,
								binBandArr);

			mChromaFilterBankArr=__vnew(num*baseNum, NULL);
			chroma_cqtFilterBank(num,baseNum,binPerOctave,
								&baseFre,
								mChromaFilterBankArr);
		}
		else if(filterScaleType==SpectralFilterBankScale_DeepChroma){ // similar cqt-chroma 
			mFilterBankArr=__vnew(num*baseNum, NULL);
			freBandArr=__vnew(baseNum+2, NULL); 
			binBandArr=__vnewi(baseNum+2, NULL);

			chroma_cqtFilterBank(num,baseNum,12,
								&baseFre,
								mFilterBankArr);

			__spectrogramObj_calDeepBandArr(spectrogramObj,freBandArr,binBandArr);
		}
		else{ // deep
			freBandArr=__vnew(num+2, NULL); 
			binBandArr=__vnewi(num+2, NULL);

			__spectrogramObj_calDeepBandArr(spectrogramObj,freBandArr,binBandArr);
		}
	}
	else{ // linear 
		freBandArr=__vnew(num+2, NULL); 
		binBandArr=__vnewi(num+2, NULL);

		__spectrogramObj_calLinearBandArr(spectrogramObj,freBandArr,binBandArr);
	}

	start=0;
	end=num-1;

	__varangei(start,end+1,1,&indexArr);
	indexLength=num;

	// 3. dct
	r=_calRadix2(num);
	if(r){ // dct加速
		fftObj_new(&fftObj, r);
	}
	else{ // 直接dct
		dctObj_new(&dctObj, num,NULL);
	}

	// 4. deep相关
	ampArr=__vnew(fftLength/2+1, NULL);
	weightArr=__vlinspace(0, samplate/2.0, fftLength/2+1, 0);
	auditory_weightA(weightArr, fftLength/2+1, NULL);

	det=samplate/(float )fftLength;
	startIndex=floorf(lowFre/det);
	endIndex=ceilf(highFre/det);

	salienceIndexArr=__vnewi(endIndex-startIndex+1, NULL);

	spectrogramObj->stftObj=stftObj;
	spectrogramObj->mFilterBankArr=mFilterBankArr;
	spectrogramObj->mChromaFilterBankArr=mChromaFilterBankArr;

	spectrogramObj->freBandArr=freBandArr;
	spectrogramObj->binBandArr=binBandArr;

	spectrogramObj->vArr1=vArr1;

	spectrogramObj->start=start;
	spectrogramObj->end=end;

	spectrogramObj->fftObj=fftObj;
	spectrogramObj->dctObj=dctObj;

	spectrogramObj->ampArr=ampArr;
	spectrogramObj->weightArr=weightArr;

	spectrogramObj->salienceIndexArr=salienceIndexArr;

	spectrogramObj->startIndex=startIndex;
	spectrogramObj->endIndex=endIndex;

	spectrogramObj->indexArr=indexArr;
	spectrogramObj->indexLength=indexLength;
}

void spectrogramObj_setDeepConfig(SpectrogramObj spectrogramObj,SpectralDeepConfig *config){
	int isRecover=0;
	int isKWeight=0;

	float maxMin=13;
	float minMax=2;

	float ratio=10;

	if(config){
		isRecover=config->isRecover;
		isKWeight=config->isKWeight;

		if(config->maxMin>config->minMax&&
			minMax>0){
			maxMin=config->maxMin;
			minMax=config->minMax;
		}

		if(config->ratio>0){
			ratio=config->ratio;
		}

		isRecover=config->isRecover;

		spectrogramObj->isRecover=isRecover;
		spectrogramObj->isKWeight=isKWeight;

		spectrogramObj->maxMin=maxMin;
		spectrogramObj->minMax=minMax;

		spectrogramObj->ratio=ratio;
	}
}

// order 1 [1,2,3,4]
void spectrogramObj_setDeepOrder(SpectrogramObj spectrogramObj,int deepOrder){

	if(deepOrder>=1&&deepOrder<=4){
		spectrogramObj->deepOrder=deepOrder;
	}
}

void spectrogramObj_setChromaDataNormalType(SpectrogramObj spectrogramObj,ChromaDataNormalType dataNormType){

	spectrogramObj->dataNormType=dataNormType;
}

void spectrogramObj_setDataNormValue(SpectrogramObj spectrogramObj,float normValue){

	if(normValue>0){
		spectrogramObj->normValue=normValue;
	}
}

int spectrogramObj_calTimeLength(SpectrogramObj spectrogramObj,int dataLength){
	int timeLength=0;

	timeLength=stftObj_calTimeLength(spectrogramObj->stftObj,dataLength);
	return timeLength;
}	


// spectrogram方法
/***
	spectrogram
	1. 涉及到mRealArr,mImageArr的缓存
	2. stft
	3. power/mag/db计算
	4. S.dot(filterBank)
****/
static void __spectrogramObj_spectrogram(SpectrogramObj spectrogramObj,float *dataArr,int dataLength,float *mRealArr1,float *mImageArr1,int nLength,int mLength,float *mDataArr,float *mPhaseArr){
	STFTObj stftObj=NULL;

	int fftLength=0; // fftLength,timeLength,num
	int timeLength=0;

	int slideLength=0;

	int num=0;
	float *mFilterBankArr=NULL;
	float *mChromaFilterBankArr=NULL;
	float *freBandArr=NULL;

	float *sumArr=NULL;
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;
	float *entropyArr=NULL;

	float *meanFreArr=NULL; 
	float *meanValueArr=NULL; 
	
	float *mRealArr=NULL;
	float *mImageArr=NULL;
	float *mSArr=NULL;

	float *energyArr=NULL;
	float *vArr1=NULL;

	int lowIndex=0; // 针对linear
	int highIndex=0;

	SpectralDataType dataType;
	SpectralFilterBankScaleType filterScaleType;

	float normValue=1;
	int baseNum=0;

	float *mCorrectFreArr=NULL; // 校正频率 timeLength*baseNum
	float *mToneFreArr=NULL;
	float *mToneFreArr1=NULL; // 临近参照

	int deepOrder=1; 

	int eFlag=0;
	int sFlag=0;
	
	stftObj=spectrogramObj->stftObj;

	filterScaleType=spectrogramObj->filterScaleType;
	dataType=spectrogramObj->dataType;

	normValue=spectrogramObj->normValue;
	baseNum=spectrogramObj->baseNum;

	fftLength=spectrogramObj->fftLength;
	
	num=spectrogramObj->num;
	mFilterBankArr=spectrogramObj->mFilterBankArr; // num*(fftLength/2+1)
	mChromaFilterBankArr=spectrogramObj->mChromaFilterBankArr; 
	freBandArr=spectrogramObj->freBandArr;

	mRealArr=spectrogramObj->mRealArr; // timeLength*fftLength
	mImageArr=spectrogramObj->mImageArr;
	mSArr=spectrogramObj->mSArr; // timeLength*(fftLength/2+1)

	energyArr=spectrogramObj->energyArr; // timeLength
	vArr1=spectrogramObj->vArr1; // fftLength

	sumArr=spectrogramObj->sumArr; // timeLength
	cArr1=spectrogramObj->cArr1;
	cArr2=spectrogramObj->cArr2;
	entropyArr=spectrogramObj->entropyArr;

	meanFreArr=spectrogramObj->meanFreArr;
	meanValueArr=spectrogramObj->meanValueArr;
	
	lowIndex=spectrogramObj->lowIndex;
	highIndex=spectrogramObj->highIndex;

	mCorrectFreArr=spectrogramObj->mCorrectFreArr;
	mToneFreArr=spectrogramObj->mToneFreArr;
	mToneFreArr1=spectrogramObj->mToneFreArr1;

	deepOrder=spectrogramObj->deepOrder;

	if(dataArr&&dataLength){
		sFlag=1;
		eFlag=1;
	}
	else{
		if(mRealArr1&&mImageArr1&&
			nLength>0&&mLength==fftLength){
			eFlag=1;
		}	
	}

	if(!eFlag){
		return;
	}

	if(sFlag){
		timeLength=stftObj_calTimeLength(stftObj,dataLength);
	}
	else{
		timeLength=nLength;
	}
	
	slideLength=spectrogramObj->slideLength;
	// det=(fftLength-slideLength)/slideLength+1; // 每次更新,避免第二次相同时溢出!!!
	if(spectrogramObj->timeLength<timeLength||
		spectrogramObj->timeLength>timeLength*2){ // 更新缓存
		free(mRealArr);
		free(mImageArr);
		free(mSArr);

		free(energyArr);

		free(sumArr);
		free(cArr1);
		free(cArr2);
		free(entropyArr);

		free(meanFreArr);
		free(meanValueArr);
		
		mRealArr=__vnew(timeLength*fftLength, NULL);
		mImageArr=__vnew(timeLength*fftLength, NULL);
		mSArr=__vnew(timeLength*(fftLength/2+1), NULL);

		energyArr=__vnew(timeLength, NULL);

		sumArr=__vnew(timeLength, NULL);
		cArr1=__vnew(timeLength, NULL);
		cArr2=__vnew(timeLength, NULL);
		entropyArr=__vnew(timeLength, NULL);

		meanFreArr=__vnew(timeLength, NULL);
		meanValueArr=__vnew(timeLength, NULL);

		if(spectrogramObj->filterScaleType==SpectralFilterBankScale_Deep||
			spectrogramObj->filterScaleType==SpectralFilterBankScale_DeepChroma){
			free(mCorrectFreArr);
			free(mToneFreArr);
			free(mToneFreArr1);

			mCorrectFreArr=__vnew(timeLength*baseNum, NULL);
			mToneFreArr=__vnew(timeLength*baseNum, NULL);
			mToneFreArr1=__vnew(timeLength*baseNum, NULL);
		}
	}

	// 针对deepFilter
	spectrogramObj->timeLength=timeLength;
	spectrogramObj->mCorrectFreArr=mCorrectFreArr;
	spectrogramObj->mToneFreArr=mToneFreArr;
	spectrogramObj->mToneFreArr1=mToneFreArr1;

	if(sFlag){
		stftObj_stft(stftObj,dataArr,dataLength,mRealArr,mImageArr);
	}
	else{
		memcpy(mRealArr, mRealArr1, sizeof(float )*nLength*fftLength);
		memcpy(mImageArr, mImageArr1, sizeof(float )*nLength*fftLength);
	}

//	// energyArr ->base stft reslut
//	for(int i=0;i<timeLength;i++){
//		__vcsquare(mRealArr+i*fftLength,mImageArr+i*fftLength,fftLength,vArr1);
//		energyArr[i]=__vsum(vArr1, fftLength)/fftLength; // logf --> mfcc
//	}

	// phase ->base stft reslut
	if(mPhaseArr&&filterScaleType==SpectralFilterBankScale_Linear){
		int len=0;

		len=highIndex-lowIndex+1;
		for(int i=0;i<timeLength;i++){
			for(int j=lowIndex,k=0;j<=highIndex;j++,k++){
				float r1=0;

				r1=mRealArr[i*fftLength+j];
				if(r1<1e-16){
					r1=1e-16;
				}

				mPhaseArr[i*len+k]=atan2f(mImageArr[i*fftLength+j], r1);
			}
		}
	}

	// 计算mR, mI => in pace mR; t*fftLength => t*(fftLength/2+1)
	if(filterScaleType==SpectralFilterBankScale_Linear&&
		spectrogramObj->lowIndex==0&&
		spectrogramObj->highIndex==spectrogramObj->fftLength/2){ // linear 0~fftLength/2 特殊处理

		__mcsquare2(mRealArr, mImageArr, timeLength, fftLength, fftLength/2+1, mDataArr); // S^2

		if(dataType==SpectralData_Mag){
			for(int i=0;i<timeLength*(fftLength/2+1);i++){
				mDataArr[i]=sqrtf(mDataArr[i]);
			}
		}
		else if(dataType==SpectralData_Power){
			if(normValue!=1){
				for(int i=0;i<timeLength*(fftLength/2+1);i++){
					mDataArr[i]=powf(mDataArr[i], normValue);
				}
			}
		}

		if(spectrogramObj->isDebug){
			printf("stft power spectrogram is :\n");
			__mdebug(mDataArr, timeLength, fftLength/2+1, 1);
			printf("\n\n");
		}
	}
	else{
		__mcsquare2(mRealArr, mImageArr, timeLength, fftLength, fftLength/2+1, mSArr); // S^2

		if(dataType==SpectralData_Mag||
			filterScaleType==SpectralFilterBankScale_Deep||
			filterScaleType==SpectralFilterBankScale_DeepChroma){ // deep/deepChroma必须先amp

			for(int i=0;i<timeLength*(fftLength/2+1);i++){
				mSArr[i]=sqrtf(mSArr[i]);
			}
		}
		else if(dataType==SpectralData_Power){
			if(normValue!=1){ 
				for(int i=0;i<timeLength*(fftLength/2+1);i++){
					mSArr[i]=powf(mSArr[i], normValue);
				}
			}
		}

		if(spectrogramObj->isDebug){
			printf("stft power spectrogram is :\n");
			__mdebug(mSArr, timeLength, fftLength/2+1, 1);
			printf("\n\n");
		}

		// S.dot(filterBank) t*(fftLength/2+1) @ num*(fftLength/2+1) => t*num
		if(filterScaleType==SpectralFilterBankScale_Mel||
			filterScaleType==SpectralFilterBankScale_Bark||
			filterScaleType==SpectralFilterBankScale_Erb||
			filterScaleType==SpectralFilterBankScale_Octave||
			filterScaleType==SpectralFilterBankScale_Linspace||
			filterScaleType==SpectralFilterBankScale_Log){ // mel/bark/erb/log linspace/logspace

			__mdot1(mSArr,mFilterBankArr,
				timeLength,fftLength/2+1,
				num,fftLength/2+1,
				mDataArr);

			if(spectrogramObj->isDebug){
				printf("filterBank is :\n");
				__mdebug(mFilterBankArr, num, fftLength/2+1, 1);
				printf("\n\n");
			}
		}
		else if(filterScaleType==SpectralFilterBankScale_Chroma){
			int type1=0;
			int p1=1;
			ChromaDataNormalType _normType=ChromaDataNormal_Max;

			if(spectrogramObj->lowIndex!=0||
				spectrogramObj->highIndex!=spectrogramObj->fftLength/2){ // 区间有效
				
				for(int i=0;i<timeLength;i++){
					for(int j=0;j<fftLength/2+1;j++){
						if(j<spectrogramObj->lowIndex||j>spectrogramObj->highIndex){
							mSArr[i*(fftLength/2+1)+j]=0;
						}
					}
				}
			}

			__mdot1(mSArr,mFilterBankArr,
				timeLength,fftLength/2+1,
				num,fftLength/2+1,
				mDataArr);

			if(dataType==SpectralData_Mag){
				if(normValue!=1){
					for(int i=0;i<timeLength*num;i++){
						mDataArr[i]=powf(mDataArr[i],normValue);
					}
				}
			}

			_normType=spectrogramObj->dataNormType;
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

				__mnormalize(mDataArr,timeLength,num,1,type1,p1,mDataArr);
			}

			if(spectrogramObj->isDebug){
				printf("filterBank is :\n");
				__mdebug(mFilterBankArr, num, fftLength/2+1, 1);
				printf("\n\n");
			}
		}
		else if(filterScaleType==SpectralFilterBankScale_LogChroma){
			int type1=0;
			int p1=1;
			ChromaDataNormalType _normType=ChromaDataNormal_Max;

			__mdot1(mSArr,mFilterBankArr,
				timeLength,fftLength/2+1,
				baseNum,fftLength/2+1,
				mImageArr);

			__mdot1(mImageArr,mChromaFilterBankArr,
				timeLength,baseNum,
				num,baseNum,
				mDataArr);

			if(dataType==SpectralData_Mag){
				if(normValue!=1){
					for(int i=0;i<timeLength*num;i++){
						mDataArr[i]=powf(mDataArr[i],normValue);
					}
				}
			}

			_normType=spectrogramObj->dataNormType;
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

				__mnormalize(mDataArr,timeLength,num,1,type1,p1,mDataArr);
			}

			if(spectrogramObj->isDebug){
				printf("filterBank is :\n");
				__mdebug(mFilterBankArr, baseNum, fftLength/2+1, 1);
				printf("\n\n");

				printf("chromaFilterBank is :\n");
				__mdebug(mFilterBankArr, num, baseNum, 1);
				printf("\n\n");
			}
		}
		else if(filterScaleType==SpectralFilterBankScale_Deep){
			int k=1;

			if(deepOrder>=1&&deepOrder<=2){
				k=3;
			}
			else{ // >=3
				k=5;
			}

			__spectrogramObj_deepFilter(spectrogramObj,mSArr,mDataArr,1);

			if(dataType==SpectralData_Power){
				for(int i=0;i<k*timeLength*num;i++){
					mDataArr[i]*=mDataArr[i];
				}

				if(normValue!=1){
					for(int i=0;i<k*timeLength*num;i++){
						mDataArr[i]=powf(mDataArr[i], normValue);
					}
				}
			}
			else if(dataType==SpectralData_Mag){
				if(normValue!=1){
					for(int i=0;i<k*timeLength*num;i++){
						mDataArr[i]=powf(mDataArr[i],normValue);
					}
				}
			}
		}
		else if(filterScaleType==SpectralFilterBankScale_DeepChroma){
			int type1=0;
			int p1=1;
			ChromaDataNormalType _normType=ChromaDataNormal_Max;

			memset(mImageArr, 0, sizeof(float )*timeLength*baseNum);
			__spectrogramObj_deepFilter(spectrogramObj,mSArr,mImageArr,0);

			if(dataType==SpectralData_Power){
				for(int i=0;i<timeLength*baseNum;i++){
					mImageArr[i]*=mImageArr[i];
				}

				if(normValue!=1){
					for(int i=0;i<timeLength*baseNum;i++){
						mImageArr[i]=powf(mImageArr[i], normValue);
					}
				}
			}

			__mdot1(mImageArr,mFilterBankArr,
				timeLength,baseNum,
				num,baseNum,
				mDataArr);

			if(dataType==SpectralData_Mag){
				if(normValue!=1){
					for(int i=0;i<timeLength*num;i++){
						mDataArr[i]=powf(mDataArr[i],normValue);
					}
				}
			}

			_normType=spectrogramObj->dataNormType;
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

				__mnormalize(mDataArr,timeLength,num,1,type1,p1,mDataArr);
			}

			if(spectrogramObj->isDebug){
				printf("filterBank is :\n");
				__mdebug(mFilterBankArr, num, fftLength/2+1, 1);
				printf("\n\n");
			}
		}
		else{ // linear情况下 
			int len=0;

			len=highIndex-lowIndex+1;
			for(int i=0;i<timeLength;i++){
				for(int j=lowIndex,k=0;j<=highIndex;j++,k++){
					mDataArr[i*len+k]=mSArr[i*(fftLength/2+1)+j];
				}
			}
		}
	}

	if(filterScaleType!=SpectralFilterBankScale_Chroma&&
		filterScaleType!=SpectralFilterBankScale_LogChroma&&
		filterScaleType!=SpectralFilterBankScale_DeepChroma&&
		filterScaleType!=SpectralFilterBankScale_Deep){ // !(chroma&deep) 

		if(dataType==SpectralData_Mag){
			if(normValue!=1){
				for(int i=0;i<timeLength*num;i++){
					mDataArr[i]=powf(mDataArr[i],normValue);
				}
			}
		}
	}

	if(spectrogramObj->isDebug){

		if(filterScaleType==SpectralFilterBankScale_Mel){
			printf("mel ");
		}
		else if(filterScaleType==SpectralFilterBankScale_Bark){
			printf("bark ");
		}
		else if(filterScaleType==SpectralFilterBankScale_Erb){
			printf("erb ");
		}
		else if(filterScaleType==SpectralFilterBankScale_Linear){
			printf("linear ");
		}
		else{
			printf("other ");
		}

		printf("spectrogram is :\n");
		__mdebug(mDataArr, timeLength, num, 1);
		printf("\n\n");
	}

	spectrogramObj->mRealArr=mRealArr;
	spectrogramObj->mImageArr=mImageArr;
	spectrogramObj->mSArr=mSArr;

	spectrogramObj->energyArr=energyArr;

	spectrogramObj->sumArr=sumArr;
	spectrogramObj->cArr1=cArr1;
	spectrogramObj->cArr2=cArr2;
	spectrogramObj->entropyArr=entropyArr;

	spectrogramObj->meanFreArr=meanFreArr;
	spectrogramObj->meanValueArr=meanValueArr;

	// 数据变换或内存变换 重置状态
	spectrogramObj->isSum=0;
	spectrogramObj->isC1=0;
	spectrogramObj->isC2=0;
	spectrogramObj->isEntropy=0;
	spectrogramObj->isEnNorm=0;

	spectrogramObj->isMean=0;
}

void spectrogramObj_spectrogram(SpectrogramObj spectrogramObj,float *dataArr,int dataLength,float *mSpectArr,float *mPhaseArr){
	
	__spectrogramObj_spectrogram(spectrogramObj,dataArr,dataLength,
								NULL,NULL,0,0,
								mSpectArr,mPhaseArr);
}

void spectrogramObj_spectrogram1(SpectrogramObj spectrogramObj,float *mRealArr,float *mImageArr,int nLength,int mLength,float *mSpectArr,float *mPhaseArr){
	
	__spectrogramObj_spectrogram(spectrogramObj,NULL,0,
								mRealArr,mImageArr,nLength,mLength,
								mSpectArr,mPhaseArr);
}

// mfcc/gtcc/???
/***
	1. log&DCT
	2. delta/deltaDelta
****/
static void __spectrogramObj_xxcc(SpectrogramObj spectrogramObj,float *mDataArr1,int mLength,CepstralRectifyType *rectifyType,float *mDataArr2){
	FFTObj fftObj=NULL;
	DCTObj dctObj=NULL;

	CepstralRectifyType recType=CepstralRectify_Log;

	int num=0;
	int timeLength=0;

	float *mRealArr=NULL; 
	float *mImageArr=NULL;

	num=spectrogramObj->num;
	timeLength=spectrogramObj->timeLength;

	mRealArr=spectrogramObj->mRealArr;
	mImageArr=spectrogramObj->mImageArr;

	fftObj=spectrogramObj->fftObj;
	dctObj=spectrogramObj->dctObj;

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

	if(spectrogramObj->isDebug){
		printf("xxcc is :\n");
		__mdebug(mDataArr2, timeLength, mLength, 1);
		printf("\n\n");
	}
}

// mel/erb ==> mfcc/gtcc
void spectrogramObj_mfcc(SpectrogramObj spectrogramObj,float *mDataArr1,int mLength,float *mDataArr2){
	SpectralFilterBankScaleType filterScaleType;

	filterScaleType=spectrogramObj->filterScaleType;
	if(filterScaleType==SpectralFilterBankScale_Mel){
		__spectrogramObj_xxcc(spectrogramObj,mDataArr1,mLength,NULL,mDataArr2);
	}
}

void spectrogramObj_gtcc(SpectrogramObj spectrogramObj,float *mDataArr1,int mLength,float *mDataArr2){
	SpectralFilterBankStyleType filterStyleType;

	filterStyleType=spectrogramObj->filterStyleType;
	if(filterStyleType==SpectralFilterBankStyle_Gammatone){
		__spectrogramObj_xxcc(spectrogramObj,mDataArr1,mLength,NULL,mDataArr2);
	}
}

void spectrogramObj_lfcc(SpectrogramObj spectrogramObj,float *mDataArr1,int mLength,float *mDataArr2){
	SpectralFilterBankScaleType filterScaleType;

	filterScaleType=spectrogramObj->filterScaleType;
	if(filterScaleType==SpectralFilterBankScale_Linear){
		__spectrogramObj_xxcc(spectrogramObj,mDataArr1,mLength,NULL,mDataArr2);
	}
}

void spectrogramObj_bfcc(SpectrogramObj spectrogramObj,float *mDataArr1,int mLength,float *mDataArr2){
	SpectralFilterBankScaleType filterScaleType;

	filterScaleType=spectrogramObj->filterScaleType;
	if(filterScaleType==SpectralFilterBankScale_Bark){
		__spectrogramObj_xxcc(spectrogramObj,mDataArr1,mLength,NULL,mDataArr2);
	}
}

void spectrogramObj_xxcc(SpectrogramObj spectrogramObj,float *mDataArr1,int mLength,CepstralRectifyType *rectifyType,float *mDataArr2){

	__spectrogramObj_xxcc(spectrogramObj,mDataArr1,mLength,rectifyType,mDataArr2);
}

/***
	mfcc standard/xxcc standard
	13*3/14*3 vector
	logEnergy?+delta+deltaDelta
	deltawindowLength 9(defalut) 必须odd>=3
	energyType Repalce
****/
void spectrogramObj_mfccStandard(SpectrogramObj spectrogramObj,float *mDataArr1,
								int *deltawindowLength,CepstralEnergyType *energyType,CepstralRectifyType *rectifyType,
								float *mDataArr2){

}

void spectrogramObj_xxccStandard(SpectrogramObj spectrogramObj,float *mDataArr1,
								int *deltawindowLength,CepstralEnergyType *energyType,CepstralRectifyType *rectifyType,
								float *mDataArr2){

}

/***
	mDataArr1 mag/power
	mDataArr2 timbre(formant) 
	mDataArr3 pitch
	只针对非chroma相关类型
****/
void spectrogramObj_deconv(SpectrogramObj spectrogramObj,float *mDataArr1,float *mDataArr2,float *mDataArr3){
	int timeLength=0;
	int num=0;
	
	FFTObj devFFTObj=NULL;
	int devFFTLength=0;

	float *devDataArr=NULL; // cqt mag&fft mag

	float *devRealArr1=NULL; // fft
	float *devImageArr1=NULL;

	float *devRealArr2=NULL; // ifft
	float *devImageArr2=NULL;

	SpectralDataType dataType=SpectralData_Mag; 

	__spectrogramObj_dealDeconv(spectrogramObj);

	dataType=spectrogramObj->dataType;

	timeLength=spectrogramObj->timeLength;
	num=spectrogramObj->num;

	// 1. 处理缓存
	__spectrogramObj_dealDeconv(spectrogramObj);
	devFFTObj=spectrogramObj->devFFTObj;
	devFFTLength=spectrogramObj->devFFTLength;

	devDataArr=spectrogramObj->devDataArr;

	devRealArr1=spectrogramObj->devRealArr1;
	devImageArr1=spectrogramObj->devImageArr1;

	devRealArr2=spectrogramObj->devRealArr2;
	devImageArr2=spectrogramObj->devImageArr2;

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

static void __spectrogramObj_dealDeconv(SpectrogramObj spectrogramObj){
	FFTObj devFFTObj=NULL;
	int devFFTLength=0;

	float *devDataArr=NULL;

	float *devRealArr1=NULL;
	float *devImageArr1=NULL;

	float *devRealArr2=NULL;
	float *devImageArr2=NULL;

	int num=0;
	int radix2Exp=0;

	num=spectrogramObj->num;
	devFFTLength=util_ceilPowerTwo(2*num);
	if(devFFTLength!=spectrogramObj->devFFTLength){
		devFFTObj=spectrogramObj->devFFTObj;

		devDataArr=spectrogramObj->devDataArr;

		devRealArr1=spectrogramObj->devRealArr1;
		devImageArr1=spectrogramObj->devImageArr1;

		devRealArr2=spectrogramObj->devRealArr2;
		devImageArr2=spectrogramObj->devImageArr2;

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

		spectrogramObj->devFFTObj=devFFTObj;
		spectrogramObj->devFFTLength=devFFTLength;

		spectrogramObj->devDataArr=devDataArr;

		spectrogramObj->devRealArr1=devRealArr1;
		spectrogramObj->devImageArr1=devImageArr1;

		spectrogramObj->devRealArr2=devRealArr2;
		spectrogramObj->devImageArr2=devImageArr2;
	}
}

/***
	deepFilter 基于amp的filter模型
	mDataArr1 timeLength*(fftLength/2+1)
	mDataArr2 timeLength*baseNum+timeLength*baseNum...
	mDataArr2 布局max+left1+right1+left2+right2

****/
static void __spectrogramObj_deepFilter(SpectrogramObj spectrogramObj,float *mDataArr1,float *mDataArr2,int isDeep){
	int fftLength=0; // fftLength,timeLength,num
	int timeLength=0;

	int samplate=0;
	int baseNum=0;

	WindowType windowType=Window_Hamm;
	int deepOrder=1; // 1

	int *salienceIndexArr=NULL; //  针对deep/deepChroma
	int salienceLength=0;

	int midiStart=0;
	int midiEnd=0;

	float det=0;

	float *mCorrectFreArr=NULL; // 校正频率 timeLength*baseNum
	float *mToneFreArr=NULL;
	float *mToneFreArr1=NULL;

	fftLength=spectrogramObj->fftLength;
	timeLength=spectrogramObj->timeLength;

	samplate=spectrogramObj->samplate;
	baseNum=spectrogramObj->baseNum;

	windowType=spectrogramObj->windowType;
	deepOrder=spectrogramObj->deepOrder;

	midiStart=spectrogramObj->midiStart;
	midiEnd=spectrogramObj->midiEnd;

	mCorrectFreArr=spectrogramObj->mCorrectFreArr;
	mToneFreArr=spectrogramObj->mToneFreArr;
	mToneFreArr1=spectrogramObj->mToneFreArr1;

	det=samplate/(float )fftLength;
	salienceIndexArr=spectrogramObj->salienceIndexArr;
	for(int i=0;i<timeLength;i++){
		float *ampArr=NULL;

		ampArr=mDataArr1+i*(fftLength/2+1);
		salienceLength=__spectrogramObj_calSalience(spectrogramObj,ampArr);
		for(int j=0;j<salienceLength;j++){
			float correctFre=0;
			float toneFre=0;
			float toneFre1=0;

			int stftIndex=0;
			int leftIndex1=0;
			int rightIndex1=0;
			int leftIndex2=0;
			int rightIndex2=0;

			int midiIndex=0; 
			int deepIndex=0;

			float cur=0;
			float left1=0;
			float right1=0;
			float left2=0;
			float right2=0;

			float scale=0;
			float value=0;

			float flag=0;

			stftIndex=salienceIndexArr[j];
			leftIndex1=stftIndex-1;
			rightIndex1=stftIndex+1;

			cur=ampArr[stftIndex];
			left1=ampArr[leftIndex1];
			right1=ampArr[rightIndex1];

			if(windowType==Window_Hamm){
				correct_hamm(cur,left1,right1,&scale,&value);
			}
			else if(windowType==Window_Hann){
				correct_hann(cur,left1,right1,&scale,&value);
			}
			else if(windowType==Window_Rect){
				correct_rect(cur,left1,right1,&scale,&value);
			}

			correctFre=(stftIndex+scale)*det;
			_calTone(correctFre, &toneFre, &toneFre1);

			midiIndex=auditory_freToMidi(toneFre);
			deepIndex=midiIndex-midiStart;

			if(deepIndex>=0&&deepIndex<baseNum){
				if(mDataArr2[i*baseNum+deepIndex]<cur){
					mDataArr2[i*baseNum+deepIndex]=cur;
					flag=1;
				}

				mCorrectFreArr[i*baseNum+deepIndex]=correctFre;
				mToneFreArr[i*baseNum+deepIndex]=toneFre;
				mToneFreArr1[i*baseNum+deepIndex]=toneFre1;
			}

			if(isDeep&&flag){ 
				if(deepOrder==1){ // left1+right1
					if(left1>right1){
						mDataArr2[timeLength*baseNum+i*baseNum+deepIndex]=left1;
					}
					else{
						mDataArr2[2*timeLength*baseNum+i*baseNum+deepIndex]=right1;
					}
				}
				else if(deepOrder==2){ // left1+right1
					mDataArr2[timeLength*baseNum+i*baseNum+deepIndex]=left1;
					mDataArr2[2*timeLength*baseNum+i*baseNum+deepIndex]=right1;
				}
				else{ // left1+right1+left2+right2
					mDataArr2[timeLength*baseNum+i*baseNum+deepIndex]=left1;
					mDataArr2[2*timeLength*baseNum+i*baseNum+deepIndex]=right1;

					leftIndex2=leftIndex1-1;
					if(leftIndex2>=0){
						left2=ampArr[leftIndex2];
						if(left2>left1){
							left2=0;
						}
					}

					rightIndex2=rightIndex1+1;
					if(rightIndex2<(fftLength/2+1)){
						right2=ampArr[rightIndex2];
						if(right2>right1){
							right2=0;
						}
					}

					if(deepOrder==3){
						if(left2>right2){
							mDataArr2[3*timeLength*baseNum+i*baseNum+deepIndex]=left2;
						}
						else{
							mDataArr2[4*timeLength*baseNum+i*baseNum+deepIndex]=right2;
						}
					}
					else{
						mDataArr2[3*timeLength*baseNum+i*baseNum+deepIndex]=left2;
						mDataArr2[4*timeLength*baseNum+i*baseNum+deepIndex]=right2;
					}
				}
			}
		}
	}
}

static int __spectrogramObj_calSalience(SpectrogramObj spectrogramObj,float *ampDataArr){
	int fftLength=0;

	int start=0; 
	int end=0;

	float maxMin=13; // 13.0
	float minMax=2; // 2.0

	float ratio=10; // 10.0

	float max=0;
	float min=0;

	int *salienceIndexArr=NULL; //  针对deep/deepChroma
	int salienceLength=0;

	start=spectrogramObj->startIndex;
	end=spectrogramObj->endIndex;

	maxMin=spectrogramObj->maxMin;
	minMax=spectrogramObj->minMax;

	ratio=spectrogramObj->ratio;
	fftLength=spectrogramObj->fftLength;

	salienceIndexArr=spectrogramObj->salienceIndexArr;
	memset(salienceIndexArr, 0, sizeof(int )*(end-start+1));

	// {
	// 	printf("ampArr is:\n");
	// 	__vdebug(ampDataArr, fftLength/2+1, 1);
	// 	printf("\n");
	// }

	__vmax(ampDataArr+start, end-start+1, &max);
	if(max<maxMin){
		return 0;
	}

	min=max/ratio;
	if(min<minMax){
		min=minMax;
	}

	if(start==0){
		start=1;
	}

	if(end==fftLength/2){
		end=fftLength/2-1;
	}

	for(int i=start;i<=end;i++){
		float _value=0;

		_value=ampDataArr[i];
		if(_value>ampDataArr[i-1]&&_value>ampDataArr[i+1]){
			if(_value>=min){
				salienceIndexArr[salienceLength]=i;
				salienceLength++;
			}
		}
	}

	spectrogramObj->salienceLength=salienceLength;
	return salienceLength;
}

// Linear/Chroma
static void __spectrogramObj_calLinearBandArr(SpectrogramObj spectrogramObj,float *freBandArr,int *binBandArr){
	float *fArr=NULL;
	int *bArr=NULL;

	int samplate=0;
	int fftLength=0;
	int lowIndex=0;

	SpectralFilterBankScaleType filterScaleType;
	int length=0;

	samplate=spectrogramObj->samplate;
	fftLength=spectrogramObj->fftLength;
	lowIndex=spectrogramObj->lowIndex;

	fArr=__vlinspace(0, samplate/2.0, fftLength/2+1, 0);
	__varangei(0, fftLength/2+1, 1, &bArr);

	filterScaleType=spectrogramObj->filterScaleType;
	if(filterScaleType==SpectralFilterBankScale_Linear){ // num
		length=spectrogramObj->num;
	}
	else if(filterScaleType==SpectralFilterBankScale_Chroma){ // lowFre~highFre
		length=spectrogramObj->baseNum;
	}

	// freBandArr+=lowIndex; ???
	memcpy(freBandArr, fArr+lowIndex, sizeof(float )*length);
	memcpy(binBandArr, bArr+lowIndex, sizeof(int )*length);

	free(fArr);
	free(bArr);
}

// Log/LogChroma
static void __spectrogramObj_calLogBandArr(SpectrogramObj spectrogramObj,float *freBandArr,int *binBandArr){
	int fftLength=0;
	int samplate=0; // filterBank 

	int binPerOctave=0; 
	float baseFre=0;

	SpectralFilterBankScaleType filterScaleType;

	float det=0;
	int index=0;
	int length=0;

	fftLength=spectrogramObj->fftLength;
	samplate=spectrogramObj->samplate;

	binPerOctave=spectrogramObj->binPerOctave;
	baseFre=spectrogramObj->baseFre;

	filterScaleType=spectrogramObj->filterScaleType;
	if(filterScaleType==SpectralFilterBankScale_Octave){ // num
		length=spectrogramObj->num;
	}
	else if(filterScaleType==SpectralFilterBankScale_LogChroma){ // lowFre~highFre
		length=spectrogramObj->baseNum;
	}

	det=samplate/(float )fftLength;
	index=auditory_freToLog(baseFre, binPerOctave);
	for(int i=index,j=0;i<index+length;i++,j++){
		float _fre=0;

		_fre=auditory_logToFre(i, binPerOctave);
		freBandArr[j]=_fre;
		binBandArr[j]=roundf(_fre/det);
	}
}

// Deep/DeepChroma
static void __spectrogramObj_calDeepBandArr(SpectrogramObj spectrogramObj,float *freBandArr,int *binBandArr){
	int fftLength=0;
	int samplate=0; // filterBank 

	float baseFre=0;

	SpectralFilterBankScaleType filterScaleType;
	
	float det=0;
	int index=0;
	int length=0;

	fftLength=spectrogramObj->fftLength;
	samplate=spectrogramObj->samplate;

	baseFre=spectrogramObj->baseFre;

	filterScaleType=spectrogramObj->filterScaleType;
	if(filterScaleType==SpectralFilterBankScale_Deep){ // num
		length=spectrogramObj->num;
	}
	else if(filterScaleType==SpectralFilterBankScale_DeepChroma){ // lowFre~highFre
		length=spectrogramObj->baseNum;
	}

	det=samplate/(float )fftLength;
	index=auditory_freToLog(baseFre, 12);
	for(int i=index,j=0;i<index+length;i++,j++){
		float _fre=0;

		_fre=auditory_logToFre(i, 12);
		freBandArr[j]=_fre;
		binBandArr[j]=roundf(_fre/det);
	}
}

// spectral相关
void spectrogramObj_setEdge(SpectrogramObj spectrogramObj,int start,int end){
	int num=0;

	num=spectrogramObj->num; // mel/bark/erb num;linear fftLength/2
	if(start>=0&&end<=num-1&&end>start){
		if(start!=spectrogramObj->start||
			end!=spectrogramObj->end){

			// 数据变换或内存变换 重置状态
			spectrogramObj->isSum=0;
			spectrogramObj->isC1=0;
			spectrogramObj->isC2=0;
			spectrogramObj->isEntropy=0;
			spectrogramObj->isEnNorm=0;

			spectrogramObj->isMean=0;
		}

		free(spectrogramObj->indexArr);

		__varangei(start, end+1, 1, &spectrogramObj->indexArr);
		spectrogramObj->indexLength=end-start+1;

		spectrogramObj->start=start;
		spectrogramObj->end=end;
	}
}

void spectrogramObj_setEdgeArr(SpectrogramObj spectrogramObj,int *indexArr,int indexLength){
	int flag=1;
	int num=0;

	num=spectrogramObj->num;
	for(int i=0;i<indexLength;i++){
		if(indexArr[i]<0||indexArr[i]>num-1){
			flag=0;
			free(indexArr);
			break;
		}
	}

	if(flag){
		spectrogramObj->isSum=0;
		spectrogramObj->isC1=0;
		spectrogramObj->isC2=0;
		spectrogramObj->isEntropy=0;
		spectrogramObj->isEnNorm=0;

		spectrogramObj->isMean=0;

		free(spectrogramObj->indexArr);

		spectrogramObj->indexArr=indexArr;
		spectrogramObj->indexLength=indexLength;

		spectrogramObj->start=indexArr[0];
		spectrogramObj->end=indexArr[indexLength-1];
	}
}

void spectrogramObj_preprocess(SpectrogramObj spectrogramObj,float *mDataArr1,float *mDataArr3){
	int num=0;
	int timeLength=0;

	float *wArr=NULL;
	int fftLength=0;

	float value=0;

	float *mArr=NULL;

	SpectralDataType dataType=SpectralData_Mag;

	num=spectrogramObj->num;
	timeLength=spectrogramObj->timeLength;

	dataType=spectrogramObj->dataType;

	if(mDataArr3){
		mArr=mDataArr3;
	}
	else {
		mArr=mDataArr1;
	}

	wArr=stftObj_getWindowDataArr(spectrogramObj->stftObj);
	fftLength=spectrogramObj->fftLength;

	value=__vsum(wArr, fftLength);
	if(dataType==SpectralData_Mag){
		value*=0.5;
	}
	else if(dataType==SpectralData_Power){
		value=0.5*value*value;
	}

	for(int i=0;i<timeLength;i++){
		for(int j=0;j<num;j++){
			mArr[i*num+j]=mDataArr1[i*num+j]/value;
			if(j==0||j==fftLength/2){
				mArr[i*num+j]*=0.5;
			}
		}
	}
}

// isSum isC1 isC2
static void __spectrogramObj_calSum(SpectrogramObj spectrogramObj,float *mDataArr){
	int num=0;
	int timeLength=0;
	
	int *indexArr=0;
	int indexLength=0;

	float *sumArr=NULL; // timeLength 内存和值都缓存

	num=spectrogramObj->num;
	timeLength=spectrogramObj->timeLength;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	sumArr=spectrogramObj->sumArr;
	for(int i=0;i<timeLength;i++){
		sumArr[i]=0; // 重置数据
		for(int j=0;j<indexLength;j++){
			int _index=0;

			_index=indexArr[j];
			sumArr[i]+=mDataArr[i*num+_index];
		}
	}

	spectrogramObj->isSum=1;
}

static void __spectrogramObj_calC1(SpectrogramObj spectrogramObj,float *mDataArr){
	int num=0;
	int timeLength=0;
	
	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *sumArr=NULL; // timeLength 内存和值都缓存
	float *cArr1=NULL; // μ1

	num=spectrogramObj->num;
	timeLength=spectrogramObj->timeLength;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;

	sumArr=spectrogramObj->sumArr;
	cArr1=spectrogramObj->cArr1;

	if(!spectrogramObj->isSum){
		__spectrogramObj_calSum(spectrogramObj,mDataArr);
	}

	spectral_centroid(mDataArr,timeLength,num,
					indexArr,indexLength,
					freArr,sumArr,
					cArr1);

	spectrogramObj->isSum=1;
	spectrogramObj->isC1=1;
}

static void __spectrogramObj_calC2(SpectrogramObj spectrogramObj,float *mDataArr){
	int num=0;
	int timeLength=0;
	
	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *sumArr=NULL; // timeLength 内存和值都缓存
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;

	num=spectrogramObj->num;
	timeLength=spectrogramObj->timeLength;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;

	sumArr=spectrogramObj->sumArr;
	cArr1=spectrogramObj->cArr1;
	cArr2=spectrogramObj->cArr2;

	if(!spectrogramObj->isSum){
		__spectrogramObj_calSum(spectrogramObj,mDataArr);
	}

	if(!spectrogramObj->isC1){
		spectral_centroid(mDataArr,timeLength,num,
					indexArr,indexLength,
					freArr,sumArr,
					cArr1);
	}

	spectral_spread(mDataArr,timeLength,num,
					indexArr,indexLength,
					freArr,sumArr,cArr1,
					cArr2);

	spectrogramObj->isSum=1;
	spectrogramObj->isC1=1;
	spectrogramObj->isC2=1;
}

static void __spectrogramObj_calEntropy(SpectrogramObj spectrogramObj,float *mDataArr,int isNorm){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *sumArr=NULL;
	float *entropyArr=NULL;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	sumArr=spectrogramObj->sumArr;
	entropyArr=spectrogramObj->entropyArr;

	if(!spectrogramObj->isSum){
		__spectrogramObj_calSum(spectrogramObj,mDataArr);
	}

	spectral_entropy(mDataArr,nLength,mLength,
					indexArr,indexLength,
					sumArr,isNorm,
					entropyArr);

	spectrogramObj->isEntropy=1;
	spectrogramObj->isEnNorm=isNorm;
}

static void __spectrogramObj_calMean(SpectrogramObj spectrogramObj,float *mDataArr){
	int num=0;
	int timeLength=0;
	
	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *meanFreArr=NULL; // 仅内存缓存
	float *meanValueArr=NULL; 

	num=spectrogramObj->num;
	timeLength=spectrogramObj->timeLength;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;

	meanFreArr=spectrogramObj->meanFreArr;
	meanValueArr=spectrogramObj->meanValueArr;
	for(int j=0;j<indexLength;j++){
		int _index=0;

		_index=indexArr[j];
		meanFreArr[0]+=freArr[_index];
	}
	meanFreArr[0]=meanFreArr[0]/indexLength;
	for(int i=1;i<timeLength;i++){
		meanFreArr[i]=meanFreArr[0];
	}

	for(int i=0;i<timeLength;i++){
		meanValueArr[i]=0;
		for(int j=0;j<indexLength;j++){
			int _index=0;

			_index=indexArr[j];
			meanValueArr[i]+=mDataArr[i*num+_index];
		}
		meanValueArr[i]=meanValueArr[i]/indexLength;
	}

	spectrogramObj->isMean=1;
}

void spectrogramObj_flatness(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;
	float *sumArr=NULL;
	
	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;
	sumArr=spectrogramObj->sumArr;

	if(!spectrogramObj->isSum){
		__spectrogramObj_calSum(spectrogramObj,mDataArr);
	}

	spectral_flatness(mDataArr,nLength,mLength,
					indexArr,indexLength,
					freArr,sumArr,
					dataArr);
}

void spectrogramObj_flux(SpectrogramObj spectrogramObj,float *mDataArr,
						int step,float p,int isPostive,int *isExp,int *type,
						float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	int _isExp=0;
	int _type=0;
	
	if(isExp){
		_isExp=*isExp;
	}

	if(type){
		_type=*type;
	}

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_flux(mDataArr,nLength,mLength,
				indexArr,indexLength,
				step,p,isPostive,_isExp,_type,
				dataArr);
}

void spectrogramObj_rolloff(SpectrogramObj spectrogramObj,float *mDataArr,float threshold,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;
	float *sumArr=NULL;
	
	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;
	sumArr=spectrogramObj->sumArr;

	if(!spectrogramObj->isSum){
		__spectrogramObj_calSum(spectrogramObj,mDataArr);
	}

	spectral_rolloff(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,sumArr,threshold,
				dataArr);
}

void spectrogramObj_centroid(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){

	if(!spectrogramObj->isC1){
		__spectrogramObj_calC1(spectrogramObj,mDataArr);
	}
	
	memcpy(dataArr, spectrogramObj->cArr1, sizeof(float )*spectrogramObj->timeLength);
}

void spectrogramObj_spread(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	
	if(!spectrogramObj->isC2){
		__spectrogramObj_calC2(spectrogramObj,mDataArr);
	}
	
	memcpy(dataArr, spectrogramObj->cArr2, sizeof(float )*spectrogramObj->timeLength);
}

void spectrogramObj_skewness(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *sumArr=NULL;
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;
	
	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;

	sumArr=spectrogramObj->sumArr;
	cArr1=spectrogramObj->cArr1;
	cArr2=spectrogramObj->cArr2;

	if(!spectrogramObj->isC2){
		__spectrogramObj_calC2(spectrogramObj,mDataArr);
	}

	spectral_skewness(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,sumArr,cArr1,cArr2,
				dataArr);
}

void spectrogramObj_kurtosis(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *sumArr=NULL;
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;
	
	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;

	sumArr=spectrogramObj->sumArr;
	cArr1=spectrogramObj->cArr1;
	cArr2=spectrogramObj->cArr2;

	if(!spectrogramObj->isC2){
		__spectrogramObj_calC2(spectrogramObj,mDataArr);
	}

	spectral_kurtosis(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,sumArr,cArr1,cArr2,
				dataArr);
}

void spectrogramObj_entropy(SpectrogramObj spectrogramObj,float *mDataArr,int isNorm,float *dataArr){

	if(!spectrogramObj->isEntropy||spectrogramObj->isEnNorm!=isNorm){
		__spectrogramObj_calEntropy(spectrogramObj, mDataArr,isNorm);
	}

	memcpy(dataArr, spectrogramObj->entropyArr, sizeof(float )*spectrogramObj->timeLength);
}

void spectrogramObj_crest(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;
	float *sumArr=NULL;
	
	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;
	sumArr=spectrogramObj->sumArr;

	if(!spectrogramObj->isSum){
		__spectrogramObj_calSum(spectrogramObj,mDataArr);
	}

	spectral_crest(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,sumArr,
				dataArr);
}

void spectrogramObj_slope(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *meanFreArr=NULL; 
	float *meanValueArr=NULL; 
	
	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;

	meanFreArr=spectrogramObj->meanFreArr;
	meanValueArr=spectrogramObj->meanValueArr;

	if(!spectrogramObj->isMean){
		__spectrogramObj_calMean(spectrogramObj,mDataArr);
	}

	spectral_slope(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,meanFreArr,meanValueArr,
				dataArr);
}

void spectrogramObj_decrease(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *sumArr=NULL;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	sumArr=spectrogramObj->sumArr;

	if(!spectrogramObj->isSum){
		__spectrogramObj_calSum(spectrogramObj,mDataArr);
	}

	spectral_decrease(mDataArr,nLength,mLength,
					indexArr,indexLength,
					sumArr,
					dataArr);
}

// p 2 >=1
void spectrogramObj_bandWidth(SpectrogramObj spectrogramObj,float *mDataArr,float p,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;
	float *cArr1=NULL; // μ1
	
	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;
	cArr1=spectrogramObj->cArr1;

	if(!spectrogramObj->isC1){
		__spectrogramObj_calC1(spectrogramObj,mDataArr);
	}

	spectral_bandWidth(mDataArr,nLength,mLength,
					indexArr,indexLength,
					freArr,cArr1,p,
					dataArr);
}

void spectrogramObj_rms(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_rms(mDataArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectrogramObj_energy(SpectrogramObj spectrogramObj,float *mDataArr,int isLog,float gamma,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	int isPower=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	if(spectrogramObj->dataType==SpectralData_Power){
		isPower=1;
	}

	spectral_energy(mDataArr,nLength,mLength,
					indexArr,indexLength,
					isPower,isLog,gamma,
					dataArr);
}

void spectrogramObj_hfc(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_hfc(mDataArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectrogramObj_sd(SpectrogramObj spectrogramObj,float *mDataArr,int step,int isPostive,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_sd(mDataArr,nLength,mLength,
				indexArr,indexLength,
				step,isPostive,
				dataArr);
}

void spectrogramObj_sf(SpectrogramObj spectrogramObj,float *mDataArr,int step,int isPostive,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_sf(mDataArr,nLength,mLength,
				indexArr,indexLength,
				step,isPostive,
				dataArr);
}

void spectrogramObj_mkl(SpectrogramObj spectrogramObj,float *mDataArr,int type,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_mkl(mDataArr,nLength,mLength,
				indexArr,indexLength,
				type,
				dataArr);
}

void spectrogramObj_pd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_pd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectrogramObj_wpd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_wpd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectrogramObj_nwpd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_nwpd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectrogramObj_cd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_cd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectrogramObj_rcd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_rcd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

// mag|power threshold >=0
void spectrogramObj_broadband(SpectrogramObj spectrogramObj,float *mDataArr,float threshold,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_broadband(mDataArr,nLength,mLength,
					indexArr,indexLength,
					threshold,
					dataArr);
}

/***
	step >=1
	threshold 0
	methodType 'sub'
	dataType 'value'
****/
void spectrogramObj_novelty(SpectrogramObj spectrogramObj,float *mDataArr,
							int step,float threshold,
							SpectralNoveltyMethodType *methodType,SpectralNoveltyDataType *dataType,
							float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	spectral_novelty(mDataArr,nLength,mLength,
					indexArr,indexLength,
					step,threshold,
					methodType,dataType,
					dataArr);
}

void spectrogramObj_eef(SpectrogramObj spectrogramObj,float *mDataArr,int isEnNorm,float *dataArr){
	int nLength=0;

	float *energyArr=NULL; 
	float *entropyArr=NULL;

	float value1=0;

	nLength=spectrogramObj->timeLength;
	
	energyArr=spectrogramObj->energyArr;
	entropyArr=spectrogramObj->entropyArr;

	if(!spectrogramObj->isEntropy||spectrogramObj->isEnNorm!=isEnNorm){
		__spectrogramObj_calEntropy(spectrogramObj,mDataArr,isEnNorm);
	}

	for(int i=0;i<nLength;i++){
		value1=energyArr[i]*entropyArr[i];
		dataArr[i]=sqrtf(1+fabsf(value1));
	}
}

void spectrogramObj_eer(SpectrogramObj spectrogramObj,float *mDataArr,int isEnNorm,float gamma,float *dataArr){
	int nLength=0;

	float *energyArr=NULL; 
	float *entropyArr=NULL;

	float value1=0;

	nLength=spectrogramObj->timeLength;

	energyArr=spectrogramObj->energyArr;
	entropyArr=spectrogramObj->entropyArr;

	if(!spectrogramObj->isEntropy||spectrogramObj->isEnNorm!=isEnNorm){
		__spectrogramObj_calEntropy(spectrogramObj,mDataArr,isEnNorm);
	}

	for(int i=0;i<nLength;i++){
		value1=logf(1+energyArr[i]*gamma)/entropyArr[i];
		dataArr[i]=sqrtf(1+fabsf(value1));
	}
}

void spectrogramObj_max(SpectrogramObj spectrogramObj,float *mDataArr,float *valueArr,float *freArr2){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float value=0;
	int index=0;
	
	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;
	for(int i=0;i<nLength;i++){
		index=indexArr[0];
		value=mDataArr[i*mLength+index];
		for(int j=1;j<indexLength;j++){
			int _index=0;

			_index=indexArr[j];
			if(value<mDataArr[i*mLength+_index]){
				index=_index;
				value=mDataArr[i*mLength+_index];
			}

		}

		valueArr[i]=value;
		freArr2[i]=freArr[index];
	}
}

void spectrogramObj_mean(SpectrogramObj spectrogramObj,float *mDataArr,float *valueArr,float *freArr2){

	if(!spectrogramObj->isMean){
		__spectrogramObj_calMean(spectrogramObj, mDataArr);
	}

	memcpy(valueArr, spectrogramObj->meanValueArr, sizeof(float )*spectrogramObj->timeLength);
	memcpy(freArr2, spectrogramObj->meanFreArr, sizeof(float )*spectrogramObj->timeLength);
}

void spectrogramObj_var(SpectrogramObj spectrogramObj,float *mDataArr,float *valueArr,float *freArr2){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *meanValueArr=NULL;
	float *meanFreArr=NULL;

	float value1=0;
	float value2=0;

	nLength=spectrogramObj->timeLength;
	mLength=spectrogramObj->num;

	indexArr=spectrogramObj->indexArr;
	indexLength=spectrogramObj->indexLength;

	freArr=spectrogramObj->freBandArr;

	meanValueArr=spectrogramObj->meanValueArr;
	meanFreArr=spectrogramObj->meanFreArr;

	if(indexLength<2){
		return;
	}

	if(!spectrogramObj->isMean){
		__spectrogramObj_calMean(spectrogramObj, mDataArr);
	}

	for(int i=0;i<nLength;i++){
		value1=0;
		value2=0;
		for(int j=0;j<indexLength;j++){
			float _value1=0;
			float _value2=0;
			int _index=0;

			_index=indexArr[j];
			_value1=meanValueArr[i]-mDataArr[i*mLength+_index];
			value1+=_value1*_value1;

			_value2=meanFreArr[i]-freArr[_index];
			value2+=_value2*_value2;
		}
		
		valueArr[i]=value1/(indexLength-1);
		freArr2[i]=value2/(indexLength-1);
	}
}

void spectrogramObj_free(SpectrogramObj spectrogramObj){
	STFTObj stftObj=NULL;

	float *mFilterBankArr=NULL;
	float *freBandArr=NULL;
	int *binBandArr=NULL;

	int *indexArr=NULL;

	float *sumArr=NULL;
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;
	float *entropyArr=NULL;

	float *meanFreArr=NULL; 
	float *meanValueArr=NULL; 
	
	float *mRealArr=NULL;
	float *mImageArr=NULL;
	float *mSArr=NULL;

	float *energyArr=NULL;
	float *vArr1=NULL;

	FFTObj fftObj=NULL;
	DCTObj dctObj=NULL;

	FFTObj devFFTObj=NULL;

	float *devDataArr=NULL; // spectral mag&fft mag

	float *devRealArr1=NULL; // fft
	float *devImageArr1=NULL;

	float *devRealArr2=NULL; // ifft
	float *devImageArr2=NULL;

	float *ampArr=NULL; // fftLength/2+1
	float *weightArr=NULL;

	int *salienceIndexArr=NULL;

	float *mCorrectFreArr=NULL; // 校正频率 timeLength*baseNum
	float *mToneFreArr=NULL;
	float *mToneFreArr1=NULL; 

	if(!spectrogramObj){
		return;
	}

	stftObj=spectrogramObj->stftObj;

	mFilterBankArr=spectrogramObj->mFilterBankArr;
	freBandArr=spectrogramObj->freBandArr;
	binBandArr=spectrogramObj->binBandArr;

	indexArr=spectrogramObj->indexArr;

	sumArr=spectrogramObj->sumArr;
	cArr1=spectrogramObj->cArr1;
	cArr2=spectrogramObj->cArr2;
	entropyArr=spectrogramObj->entropyArr;

	meanFreArr=spectrogramObj->meanFreArr;
	meanValueArr=spectrogramObj->meanValueArr;

	mRealArr=spectrogramObj->mRealArr;
	mImageArr=spectrogramObj->mImageArr;
	mSArr=spectrogramObj->mSArr;

	energyArr=spectrogramObj->energyArr;
	vArr1=spectrogramObj->vArr1;

	fftObj=spectrogramObj->fftObj;
	dctObj=spectrogramObj->dctObj;

	devFFTObj=spectrogramObj->devFFTObj;
	devDataArr=spectrogramObj->devDataArr;

	devRealArr1=spectrogramObj->devRealArr1;
	devImageArr1=spectrogramObj->devImageArr1;

	devRealArr2=spectrogramObj->devRealArr2;
	devImageArr2=spectrogramObj->devImageArr2;

	ampArr=spectrogramObj->ampArr;
	weightArr=spectrogramObj->weightArr;

	salienceIndexArr=spectrogramObj->salienceIndexArr;

	mCorrectFreArr=spectrogramObj->mCorrectFreArr;
	mToneFreArr=spectrogramObj->mToneFreArr;
	mToneFreArr1=spectrogramObj->mToneFreArr1;

	stftObj_free(stftObj);

	free(mFilterBankArr);
	free(freBandArr);
	free(binBandArr);

	free(indexArr);

	free(sumArr);
	free(cArr1);
	free(cArr2);
	free(entropyArr);

	free(meanFreArr);
	free(meanValueArr);
	
	free(mRealArr);
	free(mImageArr);
	free(mSArr);

	free(energyArr);
	free(vArr1);

	fftObj_free(fftObj);
	dctObj_free(dctObj);

	fftObj_free(devFFTObj);
	free(devDataArr);

	free(devRealArr1);
	free(devImageArr1);

	free(devRealArr2);
	free(devImageArr2);

	free(ampArr);
	free(weightArr);

	free(salienceIndexArr);

	free(mCorrectFreArr);
	free(mToneFreArr);
	free(mToneFreArr1);

	free(spectrogramObj);
}


void spectrogramObj_enableDebug(SpectrogramObj spectrogramObj,int flag){

	spectrogramObj->isDebug=1;
}

float *spectrogramObj_getFreBandArr(SpectrogramObj spectrogramObj){

	return spectrogramObj->freBandArr;
}

int *spectrogramObj_getBinBandArr(SpectrogramObj spectrogramObj){

	return spectrogramObj->binBandArr;
}

// linear => num
int spectrogramObj_getBandNum(SpectrogramObj spectrogramObj){

	return spectrogramObj->num;
}

int spectrogramObj_getBinBandLength(SpectrogramObj spectrogramObj){

	return spectrogramObj->num;
}

static void _calKWeight(float *ampArr1,int start,int end,int fftLength,float *weightArr,float *ampArr2){
	float value1=0;
	float value2=0;

	for(int i=start;i<=end;i++){
		value1=20*log10f(ampArr1[i]/fftLength);
		value1+=weightArr[i];

		value2=powf(10, value1/20)*fftLength;
		ampArr2[i]=value2;
	}
}

static int _calBaseNum(float lowFre,float highFre,float binPerOctave){
	int length=0;

	float midi1=0;
	float midi2=0;

	midi1=auditory_freToLog(lowFre,binPerOctave);
	midi2=auditory_freToLog(highFre,binPerOctave);

	length=midi2-midi1+1;

	return length;
}

static void _calTone(float value1,float *value2,float *value3){
	float _value2=0;
	float _value3=0;

	float floorIndex=0;
	float ceilIndex=0;

	float floorValue=0;
	float ceilValue=0;

	float curIndex=0;
	float nextIndex=0;
	float preIndex=0;

	float preValue=0;
	float nextValue=0;

	floorIndex=floorf(12*log2(value1/440)+69);
	ceilIndex=ceilf(12*log2(value1/440)+69);

	floorValue=powf(2,(floorIndex-69)/12)*440;
	ceilValue=powf(2,(ceilIndex-69)/12)*440;

	if(fabsf(value1-floorValue)<fabsf(value1-ceilValue)){
		curIndex=floorIndex;

		preIndex=curIndex-1;
		nextIndex=ceilIndex;

		preValue=powf(2,(preIndex-69)/12)*440;
		nextValue=powf(2,(nextIndex-69)/12)*440;
		if(fabsf(value1-preValue)<fabsf(value1-nextValue)){
			nextIndex=preIndex;
		}
	}
	else{
		curIndex=ceilIndex;
		nextIndex=floorIndex;
	}

	_value2=powf(2,(curIndex-69)/12)*440;
	_value3=powf(2,(nextIndex-69)/12)*440;

	if(value2){
		*value2=_value2;
	}

	if(value3){
		*value3=_value3;
	}
}

static float _calBaseFre(float lowFre,float binPerOctave){
	float fre=0;

	float midi=0;

	midi=auditory_freToLog(lowFre,binPerOctave);
	fre=auditory_logToFre(midi, binPerOctave);

	return fre;
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










