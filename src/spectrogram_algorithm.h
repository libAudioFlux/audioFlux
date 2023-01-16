

#ifndef SPECTROGRAM_ALGORITHM_H
#define SPECTROGRAM_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueSpectrogram *SpectrogramObj;

/***
	num>=2&&num<=2048
	samplate 32000
	lowFre linear 0 mel/bark/erb/log 27.5
	highFre 16000
	binPerOctave 12 12/24/36
	
	radix2Exp 12 
	WindowType "hann"
	slideLength 1024
	isContinue 0

	dataType "power"
	filterScaleType "linear"
	filterStyleType "slaney"
	filterNormalType "none"

	注:
		1. linear类型 num参数不起作用 有spectrogramObj_getBandNum获取 normType参数不起作用
		2. chroma类型 无spectral相关特征 normType参数不起作用 不存在freBandArr/binBandArr 
		3. mel/bark/erb minFre>=0 log/logChroma/deep/deepChroma>=32.703
		4. log/deep类型 highFre参数不起作用 
		5. binPerOctave 只针对log/logChroma类型有效
****/
int spectrogramObj_new(SpectrogramObj *spectrogramObj,int num,
					int *samplate,float *lowFre,float *highFre,int *binPerOctave,
					int *radix2Exp,WindowType *windowType,int *slideLength,int *isContinue,
					SpectralDataType *dataType,
					SpectralFilterBankScaleType *filterScaleType,
					SpectralFilterBankStyleType *filterStyleType,
					SpectralFilterBankNormalType *filterNormalType);

int spectrogramObj_newLinear(SpectrogramObj *spectrogramObj,int samplate,int radix2Exp,int *isContinue);
int spectrogramObj_newMel(SpectrogramObj *spectrogramObj,int num,int samplate,int radix2Exp,int *isContinue);
int spectrogramObj_newBark(SpectrogramObj *spectrogramObj,int num,int samplate,int radix2Exp,int *isContinue);
int spectrogramObj_newErb(SpectrogramObj *spectrogramObj,int num,int samplate,int radix2Exp,int *isContinue);
// timeLength*12 dataNormType "Max"
int spectrogramObj_newChroma(SpectrogramObj *spectrogramObj,int samplate,int radix2Exp,int *isContinue);

/***
	Deep/DeepChroma
	WindowType 'hamm' ['hamm','hann','rect']
	dataType 'power' mag||power
	deepOrder 1 [1,2,3,4]
	lowFre 32.703 highFre不起作用
	
	DeepChroma num 12
	DeepChroma dataNormType 'Max'
****/
int spectrogramObj_newDeep(SpectrogramObj *spectrogramObj,int num,int samplate,int radix2Exp,int *isContinue);
int spectrogramObj_newDeepChroma(SpectrogramObj *spectrogramObj,int samplate,int radix2Exp,int *isContinue);

// order 1 [1,2,3,4] 1/2->3*timeLength*num 3/4->5*timeLength*num
void spectrogramObj_setDeepOrder(SpectrogramObj spectrogramObj,int deepOrder);

void spectrogramObj_setChromaDataNormalType(SpectrogramObj spectrogramObj,ChromaDataNormalType dataNormType);
// 1 >0 mag=>result power=>S
void spectrogramObj_setDataNormValue(SpectrogramObj spectrogramObj,float normValue);

int spectrogramObj_calTimeLength(SpectrogramObj spectrogramObj,int dataLength);
void spectrogramObj_enableDebug(SpectrogramObj spectrogramObj,int flag);

float *spectrogramObj_getFreBandArr(SpectrogramObj spectrogramObj);
int *spectrogramObj_getBinBandArr(SpectrogramObj spectrogramObj);
// linear => num
int spectrogramObj_getBandNum(SpectrogramObj spectrogramObj);
int spectrogramObj_getBinBandLength(SpectrogramObj spectrogramObj);

// spectrogram方法 mPhaseArr can NULL; timeLength*(fftLength/2+1)
void spectrogramObj_spectrogram(SpectrogramObj spectrogramObj,float *dataArr,int dataLength,float *mSpectArr,float *mPhaseArr);
void spectrogramObj_spectrogram1(SpectrogramObj spectrogramObj,float *mRealArr,float *mImageArr,int nLength,int mLength,float *mSpectArr,float *mPhaseArr);

// mel/erb/bark ==> mfcc/gtcc/bfcc......
void spectrogramObj_mfcc(SpectrogramObj spectrogramObj,float *mDataArr1,int ccNum,float *mDataArr2);
void spectrogramObj_gtcc(SpectrogramObj spectrogramObj,float *mDataArr1,int ccNum,float *mDataArr2);
void spectrogramObj_bfcc(SpectrogramObj spectrogramObj,float *mDataArr1,int ccNum,float *mDataArr2);

/***
	mDataArr1 'mag' mag||power
	ccNum<=num
	rectifyType CepstralRectify_Log
****/
void spectrogramObj_xxcc(SpectrogramObj spectrogramObj,float *mDataArr1,int ccNum,CepstralRectifyType *rectifyType,float *mDataArr2);

/***
	mfcc standard/xxcc standard
	13*3/14*3 vector
	logEnergy?+delta+deltaDelta
	deltaWindowLength 9(defalut) 必须odd>=3
	energyType Repalce
****/
void spectrogramObj_mfccStandard(SpectrogramObj spectrogramObj,float *mDataArr1,
								int *deltaWindowLength,CepstralEnergyType *energyType,CepstralRectifyType *rectifyType,
								float *mDataArr2);
void spectrogramObj_xxccStandard(SpectrogramObj spectrogramObj,float *mDataArr1,
								int *deltaWindowLength,CepstralEnergyType *energyType,CepstralRectifyType *rectifyType,
								float *mDataArr2);

/***
	mDataArr1 'mag' mag||power
	mDataArr2 timbre(formant) 
	mDataArr3 pitch
****/
void spectrogramObj_deconv(SpectrogramObj spectrogramObj,float *mDataArr1,float *mDataArr2,float *mDataArr3);

// spectral相关 SpectralDataType=>mag&power
// start/end相对0~num-1非0~fftLength/2
void spectrogramObj_setEdge(SpectrogramObj spectrogramObj,int start,int end);
void spectrogramObj_setEdgeArr(SpectrogramObj spectrogramObj,int *indexArr,int indexLength);
void spectrogramObj_preprocess(SpectrogramObj spectrogramObj,float *mDataArr1,float *mDataArr3);

void spectrogramObj_flatness(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
// isExp 0;step 1,p 1,type 0 sum 1 mean;
void spectrogramObj_flux(SpectrogramObj spectrogramObj,float *mDataArr,
						int step,float p,int isPostive,int *isExp,int *type,
						float *dataArr);
void spectrogramObj_rolloff(SpectrogramObj spectrogramObj,float *mDataArr,float threshold,float *dataArr);

void spectrogramObj_centroid(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
void spectrogramObj_spread(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
void spectrogramObj_skewness(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
void spectrogramObj_kurtosis(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);

void spectrogramObj_entropy(SpectrogramObj spectrogramObj,float *mDataArr,int isNorm,float *dataArr);
void spectrogramObj_crest(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
void spectrogramObj_slope(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
void spectrogramObj_decrease(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
// p 2 !=0
void spectrogramObj_bandWidth(SpectrogramObj spectrogramObj,float *mDataArr,float p,float *dataArr);
void spectrogramObj_rms(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
// gamma 10
void spectrogramObj_energy(SpectrogramObj spectrogramObj,float *mDataArr,int isLog,float gamma,float *dataArr);

void spectrogramObj_hfc(SpectrogramObj spectrogramObj,float *mDataArr,float *dataArr);
void spectrogramObj_sd(SpectrogramObj spectrogramObj,float *mDataArr,int step,int isPostive,float *dataArr);
void spectrogramObj_sf(SpectrogramObj spectrogramObj,float *mDataArr,int step,int isPostive,float *dataArr);
void spectrogramObj_mkl(SpectrogramObj spectrogramObj,float *mDataArr,int type,float *dataArr);

void spectrogramObj_pd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr);
void spectrogramObj_wpd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr);
void spectrogramObj_nwpd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr);

void spectrogramObj_cd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr);
void spectrogramObj_rcd(SpectrogramObj spectrogramObj,float *mSpecArr,float *mPhaseArr,float *dataArr);

// mag|power threshold >=0
void spectrogramObj_broadband(SpectrogramObj spectrogramObj,float *mDataArr,float threshold,float *dataArr);

/***
	step >=1
	threshold 0
	methodType 'sub'
	dataType 'value'
****/
void spectrogramObj_novelty(SpectrogramObj spectrogramObj,float *mDataArr,
							int step,float threshold,
							SpectralNoveltyMethodType *methodType,SpectralNoveltyDataType *dataType,
							float *dataArr);

void spectrogramObj_eef(SpectrogramObj spectrogramObj,float *mDataArr,int isNorm,float *dataArr);
// gamma 1/10/20... song 0.5
void spectrogramObj_eer(SpectrogramObj spectrogramObj,float *mDataArr,int isNorm,float gamma,float *dataArr);

// statistics
void spectrogramObj_max(SpectrogramObj spectrogramObj,float *mDataArr,float *valueArr,float *freArr);
void spectrogramObj_mean(SpectrogramObj spectrogramObj,float *mDataArr,float *valueArr,float *freArr);
void spectrogramObj_var(SpectrogramObj spectrogramObj,float *mDataArr,float *valueArr,float *freArr);

void spectrogramObj_free(SpectrogramObj spectrogramObj);

#ifdef __cplusplus
}
#endif

#endif