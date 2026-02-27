// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../mir/_pitch_yin.h"
#include "../mir/_pitch_ffp.h"
#include "../mir/_queue.h"
#include "../mir/harmonicRatio_algorithm.h"
#include "../mir/harmonic_algorithm.h"

#include "../spectrogram_algorithm.h"

#include "tune_track.h"

struct OpaqueTuneTrack{
	int type; // 0 rule 1 queue
	int isContinue;

	PitchYINObj yinObj;
	PitchFFPObj ffpObj;
	HarmonicRatioObj hrObj;
	HarmonicObj hmObj;
	SpectrogramObj specObj;
	SpectrogramObj specObj2;

	float yinThresh; // 0.6

	float inThresh; // 0.1
	float updateThresh; // 0.5
	float cutThresh; // 0.6 

	float inFluxThresh; // 120

	int delayLength; // 1 ;1->32ms 2->64ms
	int keepLength; // 4 128ms
	int skipLength; // 3 96ms

	// status
	int inFlag;
	int keepFlag; 
	int hitFlag; // double hit
	int skipFlag;

	int inFluxFlag;
	int delayFluxLength; // 1/2

	float preFre;
	float anchorFre;
	float preValue; // minTrough
	float preDb;

	float *preFreArr;
	float *preDbArr;
	int preLength;

	float preFlux;
	float leftFlux; // preFlux2
	int onsetOffset;

	int preCount;
	int preCount2;

	float preFre2;
	float preFre3;
	float preFre4;

	// format-delay
	int formatFlag;

	float formatFre1;
	float formatFre2;
	float formatFre3;

	float formatDb1;
	float formatDb2;
	float formatDb3;

	unsigned int index;
	unsigned int formatIndex;

	int equalCount;

	// cache
	float *freArr1; 
	float *freArr2;

	float *valueArr1;
	float *valueArr2;

	float *dbArr;

	// spec
	float *mDataArr;
	float *mDataArr2;
	float *nessArr;
	float *hrArr;
	float *fluxArr;

	// spec2
	float *preSpecArr;
	int bandLength;

	// harmonic
	int *countArr;

	int timeLength;

	float updataMinValue;
	float updataMaxValue;

	int isDebug;
};

static void __tuneTrackObj_format(TuneTrackObj tuneTrackObj,float *freArr,float *dbArr,int length,float *result);

static float __updateFre2(float *freArr,float *dbArr,float *heightArr,int length,float prefre,float refFre);

static int __isSimilar(float value1,float value2);
static int __isRange(float value1,float value2);
static float __updateFre(float *arr,int length,float value,float yin,float minValue,float maxValue,int *index);
static float __compareFre(float *arr,int length,float value,int *index);

static int __isKeySimilar(float *freArr1,float *dbArr1,int length1,float *freArr2,float *dbArr2,int length2);

static float __calFlux(float *curArr,float *preArr,int length,int p,int isPostive,int isExp,int type);

static int __arr_less(float *arr1,int length,float value,float *arr2);
static int __arr_maxIndex(float *arr,int length);

/***
	samplate 32000
	lowFre 27
	highFre 4000
	radix2Exp 12
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int tuneTrackObj_new(TuneTrackObj *tuneTrackObj,
					int *samplate,float *lowFre,float *highFre,
					int *radix2Exp,int *slideLength,int *isContinue){
	int status=0;

	int _isContinue=0;

	PitchFFPObj ffpObj=NULL;
	PitchYINObj yinObj=NULL;
	HarmonicRatioObj hrObj=NULL;
	HarmonicObj hmObj=NULL;
	SpectrogramObj specObj=NULL;
	SpectrogramObj specObj2=NULL;

	SpectralDataType dataType=SpectralData_Mag; 
	SpectralFilterBankScaleType filterBankType=SpectralFilterBankScale_Linear;
	WindowType winType=Window_Hamm;

	float yinThresh=0.6; // 0.6

	float inThresh=0.25; // 0.1/0.2/ 0.17/0.18 0.21/0.22 0.3
	float updateThresh=0.5; // 0.5
	float cutThresh=0.6; // 0.6 

	float inFluxThresh=110;

	int delayLength=1; // 0->0ms 1->32ms 2->64ms
	int keepLength=4; // 128ms 
	int skipLength=4; // 128ms

	int delayFluxLength=2; 

	float *preSpecArr=NULL;
	int bandLength=0;

	TuneTrackObj tune=NULL;

	tune=*tuneTrackObj=(TuneTrackObj )calloc(1,sizeof(struct OpaqueTuneTrack ));

	if(isContinue){
		_isContinue=*isContinue;
	}

	pitchFFPObj_new(&ffpObj,
				samplate,lowFre,highFre,
				radix2Exp,slideLength,NULL,
				isContinue);

	pitchYINObj_new(&yinObj,
				samplate,lowFre,highFre,
				radix2Exp,slideLength,NULL,
				isContinue);

	harmonicRatioObj_new(&hrObj, 
						samplate, lowFre, 
						radix2Exp,&winType, slideLength);

	harmonicObj_new(&hmObj, 
					samplate, NULL,NULL,
					radix2Exp,&winType, slideLength);

	/***
		<0.2 in;
		0.2,<0.4 update;
		0.4,<0.6 keep;
		>=0.6 stop; 
		yin setThresh is not cutThresh
	****/
	pitchYINObj_setThresh(yinObj, yinThresh);

	spectrogramObj_new(&specObj,0,
					samplate,NULL,NULL,NULL,
					radix2Exp,&winType,slideLength,isContinue,
					&dataType,&filterBankType,NULL,NULL);

	{
		float _lowFre=0;
		float _highFre=400;

		spectrogramObj_new(&specObj2,0,
						samplate,&_lowFre,&_highFre,NULL,
						radix2Exp,&winType,slideLength,isContinue,
						&dataType,&filterBankType,NULL,NULL);

		bandLength=spectrogramObj_getBinBandLength(specObj2);
		preSpecArr=__vnew(bandLength, NULL);
	}


	tune->isContinue=_isContinue;

	tune->ffpObj=ffpObj;
	tune->yinObj=yinObj;
	tune->hrObj=hrObj;
	tune->hmObj=hmObj;
	tune->specObj=specObj;
	tune->specObj2=specObj2;

	tune->yinThresh=yinThresh;

	tune->inThresh=inThresh;
	tune->updateThresh=updateThresh;
	tune->cutThresh=cutThresh;

	tune->inFluxThresh=inFluxThresh;

	tune->delayLength=delayLength;
	tune->keepLength=keepLength;
	tune->skipLength=skipLength;

	tune->delayFluxLength=delayFluxLength;

	tune->preSpecArr=preSpecArr;
	tune->bandLength=bandLength;

	tune->preFreArr=__vnew(200, NULL);
	tune->preDbArr=__vnew(200, NULL);

	tune->updataMinValue=5;
	tune->updataMaxValue=8;

	return status;
}

int tuneTrackObj_calTimeLength(TuneTrackObj tuneTrackObj,int dataLength){
	int length=0;

	length=pitchFFPObj_calTimeLength(tuneTrackObj->ffpObj, dataLength);
	return length;
}

void tuneTrackObj_clear(TuneTrackObj tuneTrackObj){

	tuneTrackObj->inFlag=0;
	tuneTrackObj->preFre=0;
}

void tuneTrackObj_setTempBase(TuneTrackObj tuneTrackObj,float tempBase){

	if(tuneTrackObj->ffpObj){
		pitchFFPObj_setTempBase(tuneTrackObj->ffpObj,tempBase);
	}
}

void tuneTrackObj_setUpdateBase(TuneTrackObj tuneTrackObj,float updataMinValue,float updataMaxValue){

	tuneTrackObj->updataMinValue=(updataMinValue<1?1:updataMinValue);
	tuneTrackObj->updataMaxValue=(updataMaxValue<1?1:updataMaxValue);
}

int tuneTrackObj_getDataArr(TuneTrackObj tuneTrackObj,float **valueArr,float **dbArr,int **countArr,int **flagArr,float **fluxArr){
	int timeLen=0;

	timeLen=tuneTrackObj->timeLength;
	if(valueArr){
		*valueArr=tuneTrackObj->valueArr2;
	}

	if(dbArr){
		*dbArr=tuneTrackObj->dbArr;
	}

	if(countArr){
		*countArr=tuneTrackObj->countArr;
	}

	pitchFFPObj_getFlagData(tuneTrackObj->ffpObj,flagArr);

	if(fluxArr){
		*fluxArr=tuneTrackObj->fluxArr;
	}

	return timeLen;
}

int tuneTrackObj_getTemporalDataArr(TuneTrackObj tuneTrackObj,float **avgArr,float **maxArr,float **percentArr){
	int timeLen=0;

	timeLen=tuneTrackObj->timeLength;
	
	pitchFFPObj_getTemporalData(tuneTrackObj->ffpObj,avgArr,maxArr,percentArr);

	return timeLen;
}

void tuneTrackObj_tune(TuneTrackObj tuneTrackObj,float *dataArr,int dataLength,
					float *freArr){
	int timeLength=0;

	float *freArr1=NULL;
	float *freArr2=NULL;

	float *valueArr1=NULL;
	float *valueArr2=NULL;

	float *dbArr=NULL;

	float *mDataArr=NULL;
	float *mDataArr2=NULL;
	float *nessArr=NULL;
	float *hrArr=NULL;
	float *fluxArr=NULL;

	float *preSpecArr=NULL;
	int bandLength=0;

	int *countArr=NULL;

	float inThresh=0; // 0.2/0.1
	float updateThresh=0; // 0.4/0.5
	float cutThresh=0; // 0.6

	float inFluxThresh=0; // 120

	int mLength1=0;
	float *mFreArr=NULL;
	float *mTroughArr=NULL;
	int *lenArr1=NULL;

	int mLength2=0;
	float *mCorrArr=NULL;
	float *mDbArr=NULL;
	float *mHeightArr=NULL;	
	int *lenArr2=NULL;

	int mLength3=0;
	float *mCorrArr2=NULL;
	float *mDbArr2=NULL;
	float *mHeightArr2=NULL;	
	int *lenArr3=NULL;

	float *lightArr=NULL;

	int *formatFlagArr=NULL;

	float *formatFreArr1=NULL;
	float *formatFreArr2=NULL;
	float *formatFreArr3=NULL;

	float *formatDbArr1=NULL;
	float *formatDbArr2=NULL;
	float *formatDbArr3=NULL;

	int num=0;

	float cutFre=2000;
	float subFre=2;

	float minValue=0;
	float maxValue=0;

	preSpecArr=tuneTrackObj->preSpecArr;
	bandLength=tuneTrackObj->bandLength;

	num=spectrogramObj_getBinBandLength(tuneTrackObj->specObj);
	timeLength=pitchFFPObj_calTimeLength(tuneTrackObj->ffpObj, dataLength);

	minValue=tuneTrackObj->updataMinValue;
	maxValue=tuneTrackObj->updataMaxValue;

	if(timeLength){
		if(tuneTrackObj->timeLength<timeLength||
			tuneTrackObj->timeLength>2*timeLength){

			free(tuneTrackObj->freArr1);
			free(tuneTrackObj->freArr2);

			free(tuneTrackObj->valueArr1);
			free(tuneTrackObj->valueArr2);

			free(tuneTrackObj->dbArr);

			free(tuneTrackObj->mDataArr);
			free(tuneTrackObj->mDataArr2);
			free(tuneTrackObj->nessArr);
			free(tuneTrackObj->hrArr);
			free(tuneTrackObj->fluxArr);

			free(tuneTrackObj->countArr);

			tuneTrackObj->freArr1=__vnew(timeLength, NULL);
			tuneTrackObj->freArr2=__vnew(timeLength, NULL);

			tuneTrackObj->valueArr1=__vnew(timeLength, NULL);
			tuneTrackObj->valueArr2=__vnew(timeLength, NULL);

			tuneTrackObj->dbArr=__vnew(timeLength, NULL);

			tuneTrackObj->mDataArr=__vnew(timeLength*num, NULL);
			tuneTrackObj->mDataArr2=__vnew(timeLength*bandLength, NULL);
			tuneTrackObj->nessArr=__vnew(timeLength+1, NULL);
			tuneTrackObj->hrArr=__vnew(timeLength+1, NULL);
			tuneTrackObj->fluxArr=__vnew(timeLength, NULL);

			tuneTrackObj->countArr=__vnewi(timeLength, NULL);
		}

		freArr1=tuneTrackObj->freArr1;
		freArr2=tuneTrackObj->freArr2;

		valueArr1=tuneTrackObj->valueArr1;
		valueArr2=tuneTrackObj->valueArr2;

		dbArr=tuneTrackObj->dbArr;

		mDataArr=tuneTrackObj->mDataArr;
		mDataArr2=tuneTrackObj->mDataArr2;
		nessArr=tuneTrackObj->nessArr;
		hrArr=tuneTrackObj->hrArr;
		fluxArr=tuneTrackObj->fluxArr;

		countArr=tuneTrackObj->countArr;

		inThresh=tuneTrackObj->inThresh;
		updateThresh=tuneTrackObj->updateThresh;
		cutThresh=tuneTrackObj->cutThresh;

		inFluxThresh=tuneTrackObj->inFluxThresh;

		memset(freArr1, 0, sizeof(float )*timeLength);
		memset(freArr2, 0, sizeof(float )*timeLength);
		memset(valueArr1, 0, sizeof(float )*timeLength);
		memset(valueArr2, 0, sizeof(float )*timeLength);
		memset(dbArr, 0, sizeof(float )*timeLength);

		memset(mDataArr, 0, sizeof(float )*timeLength*num);
		memset(mDataArr, 0, sizeof(float )*timeLength*bandLength);
		memset(nessArr, 0, sizeof(float )*(timeLength+1));
		memset(hrArr, 0, sizeof(float )*(timeLength+1));
		memset(fluxArr, 0, sizeof(float )*timeLength);

		memset(countArr, 0, sizeof(int )*timeLength);

		pitchYINObj_pitch(tuneTrackObj->yinObj, dataArr, dataLength, freArr1, valueArr1, valueArr2);
		pitchFFPObj_pitch(tuneTrackObj->ffpObj, dataArr, dataLength, freArr2,dbArr);

		mLength1=pitchYINObj_getTroughData(tuneTrackObj->yinObj, &mFreArr,&mTroughArr, &lenArr1);
		mLength2=pitchFFPObj_getCorrData(tuneTrackObj->ffpObj, &mCorrArr,&mDbArr,&mHeightArr, &lenArr2);
		mLength3=pitchFFPObj_getCutData(tuneTrackObj->ffpObj, &mCorrArr2,&mDbArr2,&mHeightArr2, &lenArr3);

		pitchFFPObj_getLightData(tuneTrackObj->ffpObj,&lightArr);
		pitchFFPObj_getFormatData(tuneTrackObj->ffpObj,
								&formatFlagArr,
								&formatFreArr1,&formatFreArr2,&formatFreArr3,
								&formatDbArr1,&formatDbArr2,&formatDbArr3);

		harmonicRatioObj_harmonicRatio(tuneTrackObj->hrObj,dataArr,dataLength,hrArr);

		harmonicObj_exec(tuneTrackObj->hmObj,dataArr,dataLength);
		harmonicObj_harmonicCount(tuneTrackObj->hmObj,80,16000,countArr);

		spectrogramObj_spectrogram(tuneTrackObj->specObj,dataArr,dataLength,mDataArr,NULL);
		spectrogramObj_preprocess(tuneTrackObj->specObj,mDataArr,NULL);
		spectrogramObj_flatness(tuneTrackObj->specObj,mDataArr,nessArr);

		spectrogramObj_spectrogram(tuneTrackObj->specObj2,dataArr,dataLength,mDataArr2,NULL);

		for(int i=0;i<timeLength;i++){
			float anchorFre=0;

			tuneTrackObj->index++;

			if(tuneTrackObj->onsetOffset){
				tuneTrackObj->onsetOffset++;
			}

			fluxArr[i]=__calFlux(mDataArr2+i*bandLength,preSpecArr,bandLength,1,1,0,0);
			if(!tuneTrackObj->preFlux){
				fluxArr[i]=1e-5;
			}

			if(tuneTrackObj->inFluxFlag){
				tuneTrackObj->inFluxFlag++;
			}
			else{
				if(fluxArr[i]<tuneTrackObj->preFlux&&
					tuneTrackObj->preFlux>tuneTrackObj->leftFlux&&
					tuneTrackObj->preFlux>inFluxThresh&&
					(!tuneTrackObj->onsetOffset||tuneTrackObj->onsetOffset>5)){ // >=5*0.032

					if(fluxArr[i]>tuneTrackObj->leftFlux){
						tuneTrackObj->delayFluxLength=2;
					}
					else{
						tuneTrackObj->delayFluxLength=1;
					}

					tuneTrackObj->inFluxFlag=1;
				}
			}

			if(tuneTrackObj->inFluxFlag==tuneTrackObj->delayFluxLength){
				tuneTrackObj->inFluxFlag=0;
				tuneTrackObj->onsetOffset=1;
			}

			if(tuneTrackObj->inFlag==tuneTrackObj->delayLength+1){ // runloop

				tuneTrackObj->formatFlag=0; // reset
				tuneTrackObj->index=0; 
				tuneTrackObj->formatIndex=0;

				if(valueArr2[i]<0.2){ // <0.2 update
					if(dbArr[i]-tuneTrackObj->preDb>4&&
						!__isKeySimilar(tuneTrackObj->preFreArr, tuneTrackObj->preDbArr,tuneTrackObj->preLength,
										mCorrArr+i*mLength2, mDbArr+i*mLength2,lenArr2[i])){

						tuneTrackObj->inFlag=0;
						tuneTrackObj->keepFlag=0;
						tuneTrackObj->anchorFre=0;

						freArr[i]=tuneTrackObj->preFre;
					}
					else{
						if(tuneTrackObj->preFre<cutFre){
							freArr[i]=__updateFre(mFreArr+i*mLength1, lenArr1[i], tuneTrackObj->preFre,valueArr2[i],
												tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
							
							if(!freArr[i]&&
								tuneTrackObj->preFre>230&&tuneTrackObj->preFre<255&&
								fabsf(tuneTrackObj->preFre-freArr2[i])<15){

								freArr[i]=__updateFre(mFreArr+i*mLength1, lenArr1[i], freArr2[i],valueArr2[i],
													tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
							}
						}
						else{
							freArr[i]=__updateFre2(mCorrArr+i*mLength2,mDbArr+i*mLength2,mHeightArr+i*mLength2, lenArr2[i], tuneTrackObj->preFre,freArr2[i]);
						}

						if(freArr[i]){
							tuneTrackObj->preFre=freArr[i];
							tuneTrackObj->anchorFre=0; // reset

							tuneTrackObj->equalCount=0;
						}
						else{ // ??? hit double
							anchorFre=__updateFre(mFreArr+i*mLength1, lenArr1[i], tuneTrackObj->anchorFre,valueArr2[i],
												tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
							if(anchorFre){
								freArr[i]=tuneTrackObj->preFre;
								tuneTrackObj->anchorFre=anchorFre;
							}
							else{
								freArr[i]=tuneTrackObj->preFre;

								tuneTrackObj->equalCount++;
								// if(tuneTrackObj->equalCount>2){
								// 	tuneTrackObj->inFlag=0;
								// 	tuneTrackObj->keepFlag=0;
								// 	tuneTrackObj->anchorFre=0;
								// }
							}
						}
					}
				}
				else if(valueArr2[i]<updateThresh){ // <0.4 update ->0.2~0.4
					if(dbArr[i]-tuneTrackObj->preDb>4&&
						!__isKeySimilar(tuneTrackObj->preFreArr, tuneTrackObj->preDbArr,tuneTrackObj->preLength,
										mCorrArr+i*mLength2, mDbArr+i*mLength2,lenArr2[i])){ // change  &&!__isSimilar(tuneTrackObj->preCorr1, mCorrArr[i*mLength2])
							
						tuneTrackObj->inFlag=0;
						tuneTrackObj->keepFlag=0;
						tuneTrackObj->anchorFre=0;

						freArr[i]=tuneTrackObj->preFre;
					}
					else{
						if(tuneTrackObj->preFre<cutFre){
							freArr[i]=__updateFre(mFreArr+i*mLength1, lenArr1[i], tuneTrackObj->preFre,valueArr2[i],
												tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
						}
						else{
							freArr[i]=__updateFre2(mCorrArr+i*mLength2,mDbArr+i*mLength2,mHeightArr+i*mLength2, lenArr2[i], tuneTrackObj->preFre,freArr2[i]);
						}

						if(!freArr[i]&&valueArr2[i]>0.3){
							int _flag=0;

							_flag=__isSimilar(tuneTrackObj->preFre,freArr2[i]);
							if(_flag){
								if(fabsf(tuneTrackObj->preFre-freArr2[i])<6){
									freArr[i]=freArr2[i];
								}
							}
							else{
								_flag=__isSimilar(tuneTrackObj->preFre,freArr2[i]/2);
								if(fabsf(tuneTrackObj->preFre-freArr2[i]/2)<6){
									freArr[i]=freArr2[i]/2;
								}
							}
						}

						if(freArr[i]){
							tuneTrackObj->keepFlag=0;
							tuneTrackObj->preFre=freArr[i];
							tuneTrackObj->anchorFre=0; // reset

							tuneTrackObj->equalCount=0;
						}
						else{ // keep 4 128ms
							anchorFre=__updateFre(mFreArr+i*mLength1, lenArr1[i], tuneTrackObj->anchorFre,valueArr2[i],
												tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
							if(anchorFre){
								freArr[i]=tuneTrackObj->preFre;
								tuneTrackObj->anchorFre=anchorFre;
							}
							else{
								freArr[i]=tuneTrackObj->preFre;

								tuneTrackObj->keepFlag++;
								if(tuneTrackObj->keepFlag>tuneTrackObj->keepLength){
									tuneTrackObj->inFlag=0;
									tuneTrackObj->keepFlag=0;
									tuneTrackObj->anchorFre=0;
								}

								tuneTrackObj->equalCount++;
								// if(tuneTrackObj->equalCount>2){
								// 	tuneTrackObj->inFlag=0;
								// 	tuneTrackObj->keepFlag=0;
								// 	tuneTrackObj->anchorFre=0;
								// }
							}
						}
					}
				}
				else if(valueArr2[i]<cutThresh){ // <0.6 keep ->0.4~0.6
					if(dbArr[i]-tuneTrackObj->preDb>4){ // change  
							
						tuneTrackObj->inFlag=0;
						tuneTrackObj->keepFlag=0;
						tuneTrackObj->anchorFre=0;

						freArr[i]=tuneTrackObj->preFre;
					}
					else{
						if(tuneTrackObj->preFre<cutFre){
							freArr[i]=__updateFre(mFreArr+i*mLength1, lenArr1[i], tuneTrackObj->preFre,valueArr2[i],
												tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
						}
						else{
							freArr[i]=__updateFre2(mCorrArr+i*mLength2,mDbArr+i*mLength2,mHeightArr+i*mLength2, lenArr2[i], tuneTrackObj->preFre,freArr2[i]);
						}

						if(!freArr[i]){
							int _flag=0;

							_flag=__isSimilar(tuneTrackObj->preFre,freArr2[i]);
							if(_flag){
								if(fabsf(tuneTrackObj->preFre-freArr2[i])<6){
									freArr[i]=freArr2[i];
								}
							}
							else{
								_flag=__isSimilar(tuneTrackObj->preFre,freArr2[i]/2);
								if(fabsf(tuneTrackObj->preFre-freArr2[i]/2)<6){
									freArr[i]=freArr2[i]/2;
								}
							}
						}

						if(freArr[i]){
							tuneTrackObj->keepFlag=0;

							if(tuneTrackObj->preFre<cutFre){
								freArr[i]=tuneTrackObj->preFre;
							}
							else{
								tuneTrackObj->preFre=freArr[i];
							}

							tuneTrackObj->anchorFre=0; // reset

							tuneTrackObj->equalCount=0;
						}
						else{ // keep 4 128ms
							anchorFre=__updateFre(mFreArr+i*mLength1, lenArr1[i], tuneTrackObj->anchorFre,valueArr2[i],
												tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
							if(anchorFre){
								freArr[i]=tuneTrackObj->preFre;
								tuneTrackObj->anchorFre=anchorFre;
							}
							else{
								freArr[i]=tuneTrackObj->preFre;

								tuneTrackObj->keepFlag++;
								if(tuneTrackObj->keepFlag>tuneTrackObj->keepLength){
									tuneTrackObj->inFlag=0;
									tuneTrackObj->keepFlag=0;
									tuneTrackObj->anchorFre=0;
								}

								tuneTrackObj->equalCount++;
								// if(tuneTrackObj->equalCount>2){
								// 	tuneTrackObj->inFlag=0;
								// 	tuneTrackObj->keepFlag=0;
								// 	tuneTrackObj->anchorFre=0;
								// }
							}
						}
					}
				}
				else{ // stop 

					if(dbArr[i]-tuneTrackObj->preDb>4){ // change  
							
						tuneTrackObj->inFlag=0;
						tuneTrackObj->keepFlag=0;
						tuneTrackObj->anchorFre=0;

						freArr[i]=tuneTrackObj->preFre;
					}
					else{
						if(tuneTrackObj->preFre<cutFre){
							freArr[i]=__updateFre(mFreArr+i*mLength1, lenArr1[i], tuneTrackObj->preFre,valueArr2[i],
												tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
						}
						else{
							freArr[i]=__updateFre2(mCorrArr+i*mLength2,mDbArr+i*mLength2,mHeightArr+i*mLength2, lenArr2[i], tuneTrackObj->preFre,freArr2[i]);
						}

						if(freArr[i]){
							tuneTrackObj->keepFlag=0;

							if(tuneTrackObj->preFre<cutFre){
								freArr[i]=tuneTrackObj->preFre;
							}
							else{
								tuneTrackObj->preFre=freArr[i];
							}

							tuneTrackObj->anchorFre=0; // reset

							tuneTrackObj->equalCount=0;
						}
						else{ // keep 4 128ms
							anchorFre=__updateFre(mFreArr+i*mLength1, lenArr1[i], tuneTrackObj->anchorFre,valueArr2[i],
												tuneTrackObj->updataMinValue, tuneTrackObj->updataMaxValue,NULL);
							if(anchorFre){
								freArr[i]=tuneTrackObj->preFre;
								tuneTrackObj->anchorFre=anchorFre;
							}
							else{
								freArr[i]=tuneTrackObj->preFre;

								tuneTrackObj->keepFlag++;
								if(tuneTrackObj->keepFlag>tuneTrackObj->keepLength){
									tuneTrackObj->inFlag=0;
									tuneTrackObj->keepFlag=0;
									tuneTrackObj->anchorFre=0;
								}

								tuneTrackObj->equalCount++;
								// if(tuneTrackObj->equalCount>2){
								// 	tuneTrackObj->inFlag=0;
								// 	tuneTrackObj->keepFlag=0;
								// 	tuneTrackObj->anchorFre=0;
								// }
							}
						}
					}
				}

				// if(tuneTrackObj->preValue>0.4&&
				// 	valueArr2[i]<inThresh){ // overflow ->stop

				// 	tuneTrackObj->inFlag=0;
				// 	tuneTrackObj->keepFlag=0;
				// 	tuneTrackObj->anchorFre=0;

				// 	freArr[i]=0;
				// }
			}
			else{ // entry
				tuneTrackObj->equalCount=0;

				if(valueArr2[i]<inThresh&&
					((valueArr2[i]<0.1&&(lightArr[i]>0.98?countArr[i]>=3:1))|| // &&(lightArr[i]>0.98?countArr[i]>=4:1)
					(valueArr2[i]>=0.1&&valueArr2[i]<0.2&&
						countArr[i]>=6&&(nessArr[i]<0.13||hrArr[i]>0.8))||
					(valueArr2[i]>=0.2&&tuneTrackObj->preValue<0.2&&
						countArr[i]>=6&&(nessArr[i]<0.12||hrArr[i]>0.8)))&&
					freArr2[i]){ // &&nessArr[i]<0.13 countArr[i]>=6/9  ->ness bug!!! 

					tuneTrackObj->inFlag++;

					if(freArr2[i]>215&&freArr2[i]<225&&
						valueArr2[i]<0.1&&
						countArr[i]<=12){ // 105~115

						int flag=0;

						for(int j=0;j<lenArr3[i];j++){
							if(mCorrArr2[i*mLength3+j]>104.5&&
								mCorrArr2[i*mLength3+j]<115){

								flag=1;
								break;
							}
							else if(mCorrArr2[i*mLength3+j]>115){
								break;
							}
						}

						if(!flag){
							for(int j=0;j<lenArr2[i];j++){
								if(mCorrArr[i*mLength2+j]>104.5&&
									mCorrArr[i*mLength2+j]<115){

									flag=1;
									break;
								}
								else if(mCorrArr[i*mLength2+j]>115){
									break;
								}
							}
						}

						if(flag||1){
							if((tuneTrackObj->preFre4>105&&tuneTrackObj->preFre4<115)||
								(tuneTrackObj->preFre3>105&&tuneTrackObj->preFre3<115)){

								freArr2[i]/=2;
							}
						}
					}

					if(tuneTrackObj->preFre4>105&&tuneTrackObj->preFre4<115){ // 110-147
						if(valueArr2[i]<0.1&&
							lenArr2[i]>10){

							float *freArr=NULL;
							float *dbArr=NULL;

							int len=0;

							int index1=0;

							int k1=0;
							int k2=0;

							freArr=mCorrArr2+i*mLength3;
							dbArr=mDbArr2+i*mLength3;

							len=lenArr3[i];

							if(freArr[0]>105&&freArr[0]<115&&
								freArr[1]>140&&freArr[1]<155&&
								dbArr[1]>dbArr[2]&&
								dbArr[1]>dbArr[3]){

								freArr2[i]=freArr[1];
							}
						}
					}
					else if(tuneTrackObj->preFre4>140&&tuneTrackObj->preFre4<155){ // 147-196
						if(((freArr2[i]>95&&freArr2[i]<103)||
							(freArr2[i]>45&&freArr2[i]<50))&&
							valueArr2[i]<0.2&&
							lenArr2[i]>10){ // lenArr2[i] can more countArr[i]!!!

							float *freArr=NULL;
							float *dbArr=NULL;

							int len=0;

							int index1=0;

							freArr=mCorrArr2+i*mLength3;
							dbArr=mDbArr2+i*mLength3;

							len=lenArr3[i];

							index1=__arr_maxIndex(dbArr, len);
							if(index1==1&&
								freArr[1]>190&&freArr[1]<205&&
								dbArr[1]-dbArr[0]>8&&
								dbArr[1]-dbArr[2]>8){

								freArr2[i]=freArr[1];
							}
							else if(index1==2&&
								freArr[2]>190&&freArr[2]<205&&
								dbArr[2]-dbArr[1]>8&&
								dbArr[2]-dbArr[3]>8){

								freArr2[i]=freArr[2];
							}
						}
					}
					else if(tuneTrackObj->preFre4>240&&tuneTrackObj->preFre4<255){ // 247-329
						if(valueArr2[i]<0.1&&
							lenArr2[i]>10){

							float *freArr=NULL;
							float *dbArr=NULL;

							int len=0;

							int index1=0;

							int k1=0;
							int k2=0;

							freArr=mCorrArr2+i*mLength3;
							dbArr=mDbArr2+i*mLength3;

							len=lenArr3[i];

							index1=__arr_maxIndex(dbArr, len);
							if(freArr[index1]>300&&freArr[index1]<360){

								freArr2[i]=freArr[index1];
							}
						}
					}

					if(freArr2[i]>50&&freArr2[i]<60&&
						valueArr2[i]>0.1){

						tuneTrackObj->inFlag--;
					}
					else if(freArr2[i]>40&&freArr2[i]<50&&
						valueArr2[i]>0.1){

						// int _index=0;

						// _index=__arr_maxIndex(mDbArr+i*mLength2, lenArr2[i]);
						// if(mCorrArr[i*mLength2+_index]>220&&
						// 	mCorrArr[i*mLength2+_index]<300){

						// 	tuneTrackObj->inFlag--;
						// }
					}
					else if(freArr2[i]>160&&freArr2[i]<170&&
						valueArr2[i]<0.1&&
						countArr[i]<=3){ // 80~85

						tuneTrackObj->inFlag--;
					}
					else if(freArr2[i]>235&&freArr2[i]<260&&
						valueArr2[i]<0.1&&
						countArr[i]<=4){ // 80~85

						if((tuneTrackObj->preFre4>75&&tuneTrackObj->preFre4<90)||
							(tuneTrackObj->preFre3>75&&tuneTrackObj->preFre3<90)){

							tuneTrackObj->inFlag=0;
						}
					}
					else if(freArr2[i]>430&&freArr2[i]<450&&
						valueArr2[i]<0.1&&
						countArr[i]<=4){ // 140~155

						if((tuneTrackObj->preFre4>140&&tuneTrackObj->preFre4<155)||
							(tuneTrackObj->preFre3>140&&tuneTrackObj->preFre3<155)){

							tuneTrackObj->inFlag=0;
						}
					}
					else if(freArr2[i]>210&&freArr2[i]<230&&
						valueArr2[i]<0.1&&
						countArr[i]<=6){ // 105~115

						if((tuneTrackObj->preFre4>105&&tuneTrackObj->preFre4<115)||
							(tuneTrackObj->preFre3>105&&tuneTrackObj->preFre3<115)){

							tuneTrackObj->inFlag=0;
						}
					}
					else if(tuneTrackObj->preFre4>240&&tuneTrackObj->preFre4<255){ // 247

						int k1=0;
						int k2=0;

						float fre1=0;
						float fre2=0;

						fre1=(tuneTrackObj->preFre4>freArr2[i]?freArr2[i]:tuneTrackObj->preFre4);
						fre2=(tuneTrackObj->preFre4>freArr2[i]?tuneTrackObj->preFre4:freArr2[i]);

						__queue_fre2(fre1, fre2, &k1, &k2);
						if(k1==1&&k2==2&&
							fabsf(fre1*2-fre2)<4){

							tuneTrackObj->inFlag=0;
						}
					}
					else if(tuneTrackObj->preFre4>320&&tuneTrackObj->preFre4<345){ // 330
						if(freArr2[i]>105&&freArr2[i]<115&&
							mCorrArr2[i*mLength3]/2>105&&mCorrArr2[i*mLength3]/2<115&&
							mHeightArr2[i*mLength3]<12&&
							lenArr2[i]<=4){

							tuneTrackObj->inFlag=0;
						}
					}

					if(freArr2[i]>230){
						subFre=5;
					}
					else{
						subFre=2;
					}

					// format
					// if(tuneTrackObj->inFlag==tuneTrackObj->delayLength+1){

					// 	if(tuneTrackObj->formatFlag){
					// 		if(tuneTrackObj->index-tuneTrackObj->formatIndex>=2&&
					// 			tuneTrackObj->index-tuneTrackObj->formatIndex<=4){

					// 			__tuneTrackObj_format(tuneTrackObj, mCorrArr2+mLength3*i,mDbArr2+mLength3*i ,lenArr3[i], freArr2+i);

					// 			// tuneTrackObj->formatFlag=0;
					// 			// tuneTrackObj->formatIndex=0;
					// 		}
					// 		else if(tuneTrackObj->index-tuneTrackObj->formatIndex>4){
					// 			tuneTrackObj->formatFlag=0;
					// 			tuneTrackObj->formatIndex=0;
					// 		}
					// 	}

					// 	if(!tuneTrackObj->formatFlag&&formatFlagArr[i]){

					// 		tuneTrackObj->formatFlag=1;
					// 		tuneTrackObj->formatIndex=tuneTrackObj->index;

					// 		tuneTrackObj->formatFre1=formatFreArr1[i];
					// 		tuneTrackObj->formatFre2=formatFreArr2[i];

					// 		tuneTrackObj->formatDb1=formatDbArr1[i];
					// 		tuneTrackObj->formatDb2=formatDbArr2[i];

					// 		tuneTrackObj->inFlag=0; // reset
					// 	}
					// }

					if(tuneTrackObj->inFlag==tuneTrackObj->delayLength+1){
						int _index=-1;
						
						freArr[i]=__compareFre(mFreArr+i*mLength1, lenArr1[i], freArr2[i],&_index);
						if(freArr[i]){
							if(fabsf(freArr2[i]-mFreArr[i*mLength1+_index])<subFre){
								freArr[i]=freArr2[i];
								tuneTrackObj->preFre=freArr[i];
							}
							else{
								freArr[i]=0;
								tuneTrackObj->inFlag--;
							}
						}
						else{
							if(lenArr1[i]&&freArr2[i]){
								if(freArr2[i]>mFreArr[i*mLength1]){ // get freArr2[i]
									freArr[i]=freArr2[i];
									tuneTrackObj->preFre=freArr[i];

									tuneTrackObj->anchorFre=mFreArr[i*mLength1];
								}
							}

							if(!tuneTrackObj->anchorFre){
								tuneTrackObj->inFlag--;
							}
						}
					}
				}
				else if(valueArr2[i]<0.16&&valueArr2[i]>0.09&&
						countArr[i]>=4&&
						lightArr[i]>0.98){

					tuneTrackObj->inFlag++;

					if(freArr2[i]>230){
						subFre=2;
					}
					else{
						subFre=2;
					}

					if(tuneTrackObj->inFlag==tuneTrackObj->delayLength+1&&
						freArr2[i]){

						int _index=-1;
						
						freArr[i]=__compareFre(mFreArr+i*mLength1, lenArr1[i], freArr2[i],&_index);
						if(freArr[i]){
							if(fabsf(freArr2[i]-mFreArr[i*mLength1+_index])<subFre){
								freArr[i]=freArr2[i];
								tuneTrackObj->preFre=freArr[i];
							}
							else{
								freArr[i]=0;
								tuneTrackObj->inFlag--;
							}
						}
						else{
							if(lenArr1[i]&&freArr2[i]){
								if(freArr2[i]>mFreArr[i*mLength1]){ // get freArr2[i]
									freArr[i]=freArr2[i];
									tuneTrackObj->preFre=freArr[i];

									tuneTrackObj->anchorFre=mFreArr[i*mLength1];
								}
							}

							if(!tuneTrackObj->anchorFre){
								tuneTrackObj->inFlag--;
							}
						}
					}
				}
				else if(valueArr2[i]<0.4&&
						(countArr[i]>9||
							(tuneTrackObj->preCount>9&&
								tuneTrackObj->preCount2>9))&&
						lightArr[i]>0.98){

					tuneTrackObj->inFlag++;
					tuneTrackObj->delayLength=2;

					if(freArr2[i]>230){
						subFre=5;
					}
					else{
						subFre=2;
					}

					if(tuneTrackObj->inFlag==tuneTrackObj->delayLength+1){
						if(freArr2[i]){
							int _index=-1;
						
							freArr[i]=__compareFre(mFreArr+i*mLength1, lenArr1[i], freArr2[i],&_index);
							if(freArr[i]){
								if(fabsf(freArr2[i]-mFreArr[i*mLength1+_index])<subFre){
									freArr[i]=freArr2[i];
									tuneTrackObj->preFre=freArr[i];
								}
								else{
									freArr[i]=0;
								}
							}
							else{
								if(lenArr1[i]&&freArr2[i]){
									if(freArr2[i]>mFreArr[i*mLength1]){ // get freArr2[i]
										freArr[i]=freArr2[i];
										tuneTrackObj->preFre=freArr[i];

										tuneTrackObj->anchorFre=mFreArr[i*mLength1];
									}
								}
							}

							if(freArr[i]){
								tuneTrackObj->delayLength=1;
								tuneTrackObj->inFlag=tuneTrackObj->delayLength+1;
							}
							else{
								tuneTrackObj->inFlag--;
							}
						}
						else{
							tuneTrackObj->inFlag--;
						}
					}
				}
				// else if(countArr[i]>9&&
				// 		tuneTrackObj->onsetOffset>=3&&
				// 		tuneTrackObj->onsetOffset<=4){

				// 	int _index=-1;
						
				// 	freArr[i]=__compareFre(mFreArr+i*mLength1, lenArr1[i], freArr2[i],&_index);
				// 	if(freArr[i]){
				// 		if(fabsf(freArr2[i]-mFreArr[i*mLength1+_index])<2){
				// 			freArr[i]=freArr2[i];
				// 			tuneTrackObj->preFre=freArr[i];
				// 		}
				// 		else{
				// 			freArr[i]=0;
				// 		}
				// 	}
				// 	else{
				// 		if(lenArr1[i]&&freArr2[i]){
				// 			if(freArr2[i]>mFreArr[i*mLength1]){ // get freArr2[i]
				// 				freArr[i]=freArr2[i];
				// 				tuneTrackObj->preFre=freArr[i];

				// 				tuneTrackObj->anchorFre=mFreArr[i*mLength1];
				// 			}
				// 		}
				// 	}

				// 	if(freArr[i]){
				// 		tuneTrackObj->delayLength=1;
				// 		tuneTrackObj->inFlag=tuneTrackObj->delayLength+1;
				// 	}
				// }
				else{
					tuneTrackObj->inFlag=0;
					tuneTrackObj->keepFlag=0;
					tuneTrackObj->anchorFre=0;

					tuneTrackObj->delayLength=1;
				}
			}

			tuneTrackObj->preDb=dbArr[i];
			tuneTrackObj->preValue=valueArr2[i];

			memcpy(tuneTrackObj->preFreArr, mCorrArr+i*mLength2, sizeof(float )*lenArr2[i]);
			memcpy(tuneTrackObj->preDbArr, mDbArr+i*mLength2, sizeof(float )*lenArr2[i]);
			tuneTrackObj->preLength=lenArr2[i];

			memcpy(preSpecArr, mDataArr2+i*bandLength, sizeof(float )*bandLength);
			tuneTrackObj->leftFlux=tuneTrackObj->preFlux;
			tuneTrackObj->preFlux=fluxArr[i];

			tuneTrackObj->preCount2=tuneTrackObj->preCount;
			tuneTrackObj->preCount=countArr[i];

			tuneTrackObj->preFre4=tuneTrackObj->preFre3;
			tuneTrackObj->preFre3=tuneTrackObj->preFre2;
			tuneTrackObj->preFre2=tuneTrackObj->preFre;
		}
	}
	
	tuneTrackObj->timeLength=timeLength;
}

static void __tuneTrackObj_format(TuneTrackObj tuneTrackObj,float *freArr,float *dbArr,int length,float *result){
	float fre1=0;
	float fre2=0;

	float db1=0;
	float db2=0;

	int flag1=0;
	int flag2=0;

	int index1=0;
	int index2=0;

	int k1=0;
	int k2=0;

	float det1=0;
	float det2=0;

	fre1=tuneTrackObj->formatFre1;
	fre2=tuneTrackObj->formatFre2;

	db1=tuneTrackObj->formatDb1;
	db2=tuneTrackObj->formatDb2;

	for(int i=0;i<length;i++){
		if(fabsf(fre1-freArr[i])<10){
			flag1=1;
			index1=i;
			break;
		}
	}

	for(int i=0;i<length;i++){
		if(fabsf(fre2-freArr[i])<10){
			flag2=1;
			index2=i;
			break;
		}
	}

	if(flag1&&!flag2){
		*result=freArr[index1];
	}

	if(flag1&&flag2){
		__queue_fre2(freArr[index1], freArr[index2], &k1, &k2);
		if(k1==1&&(k2==2||k2==3)){
			det1=dbArr[index1]-db1;
			det2=dbArr[index2]-db2;

			if(det1>0&&det1>det2){
				*result=freArr[index2]/k2;
			}
			else if(det1<0&&det1>det2){
				*result=freArr[index2]/k2;
			}
		}
	}
}

void tuneTrackObj_enableDebug(TuneTrackObj tuneTrackObj,int isDebug){

	tuneTrackObj->isDebug=isDebug;
}

void tuneTrackObj_free(TuneTrackObj tuneTrackObj){

	if(tuneTrackObj){
		pitchFFPObj_free(tuneTrackObj->ffpObj);
		pitchYINObj_free(tuneTrackObj->yinObj);
		harmonicRatioObj_free(tuneTrackObj->hrObj);
		harmonicObj_free(tuneTrackObj->hmObj);
		spectrogramObj_free(tuneTrackObj->specObj);
		spectrogramObj_free(tuneTrackObj->specObj2);

		free(tuneTrackObj->freArr1);
		free(tuneTrackObj->freArr2);

		free(tuneTrackObj->valueArr1);
		free(tuneTrackObj->valueArr2);

		free(tuneTrackObj->dbArr);

		free(tuneTrackObj->mDataArr);
		free(tuneTrackObj->mDataArr2);
		free(tuneTrackObj->nessArr);
		free(tuneTrackObj->hrArr);
		free(tuneTrackObj->fluxArr);

		free(tuneTrackObj->preSpecArr);

		free(tuneTrackObj->countArr);

		free(tuneTrackObj->preFreArr);
		free(tuneTrackObj->preDbArr);

		free(tuneTrackObj);
	}
}

static int __isSimilar(float value1,float value2){
	int flag=0;

	int midi1=0;
	int midi2=0;

	midi1=util_freToMidi(value1);
	midi2=util_freToMidi(value2);
	if(abs(midi1-midi2)<=1){
		flag=1;
	}

	return flag;
}

static int __isRange(float value1,float value2){
	int flag=0;

	int midi1=0;
	int midi2=0;

	midi1=util_freToMidi(value1);
	midi2=util_freToMidi(value2);
	if(abs(midi1-midi2)<=2){
		flag=1;
	}

	return flag;
}

static float __updateFre2(float *freArr,float *dbArr,float *heightArr,int length,float preFre,float refFre){
	float fre=0;

	float value=0;

	int _index=0;

	if(!length){
		return 0;
	}

	// 1. direct
	value=fabsf(preFre-refFre);
	if(value<10){ // <5???
		return refFre;
	}

	// 2. search
	for(int i=0;i<length;i++){
		value=fabsf(freArr[i]-preFre);
		if(value<10){
			return freArr[i];
		}
	}

	// 3. times
	_index=__arr_maxIndex(dbArr, length);
	if(heightArr[_index]>15){
		for(int i=2;i<10;i++){
			value=fabsf(freArr[_index]/i-preFre);
			if(value<10){
				return freArr[_index]/i;
			}
		}
	}

	return fre;
}

static float __updateFre(float *arr,int length,float value,float yin,float minValue,float maxValue,int *index){
	float fre=0;
	int _index=-1;

	float error=5000;

	float sub=5; // 5
	float sub2=10;

	int flag=0;

	sub=minValue;
	if(value>220){
		sub=8; // 8
		sub=maxValue;
	}

	if(yin>0.3){
		sub2=minValue;
	}

	if(!length){
		return 0;
	}

	for(int i=0;i<length;i++){
		float _value=0;

		_value=fabsf(arr[i]-value);
		// if(_value<sub){ // 5/6/8
		// 	if(error>_value){
		// 		error=_value;
				
		// 		fre=arr[i];
		// 		_index=i;
		// 	}
		// }

		if(error>_value){
			error=_value;
			
			fre=arr[i];
			_index=i;
		}
	}

	if(arr[_index]>value){
		if(error<sub){
			flag=1;
		}
	}
	else{
		if(error<sub2){
			flag=1;
		}
	}

	if(!flag){
		fre=0;
	}

	if(index){
		*index=_index;
	}
	
	return fre;
}

static float __compareFre(float *arr,int length,float value,int *index){
	float fre=0;
	int _index=-1;

	float error=100;

	for(int i=0;i<length;i++){
		float _value=0;

		_value=fabsf(arr[i]-value);
		if(__isSimilar(arr[i], value)){
			if(error>_value){
				error=_value;

				fre=arr[i];
				_index=i;
			}
		}
	}

	if(index){
		*index=_index;
	}

	return fre;
}

// only for twist decay
static int __isKeySimilar(float *freArr1,float *dbArr1,int length1,float *freArr2,float *dbArr2,int length2){
	int flag=0;

	float db1=0;
	float db2=0;

	if(length1>1&&length2>1&&
		length2<=6){
		// dB desc
		__vcorrsort1(dbArr1,freArr1,NULL,NULL,length1,1);
		__vcorrsort1(dbArr2,freArr2,NULL,NULL,length2,1);

		db1=dbArr1[0];
		db2=dbArr2[0];
		if(fabsf(db1-db2)>5.6){
			return 0;
		}

		// fre asc
		__vcorrsort1(freArr1,dbArr1,NULL,NULL,2,0);
		__vcorrsort1(freArr2,dbArr2,NULL,NULL,2,0);

		flag=1;
		for(int i=0;i<2;i++){
			int k=0;

			k=util_calRangeTimes(freArr1[i], freArr2[i], NULL);
			if(k!=1){
				flag=0;
				break;
			}
		}

		if(!flag&&
			length2==2&&length1<=3){ // {82,164} ->{89,190}

			int k1=0;
			int k2=0;

			__queue_fre2(freArr1[0],freArr1[1],
						&k1,&k2);

			if(k1==1&&k2==2&&
				fabsf(freArr1[0]*2-freArr1[1])<5){

				if(freArr2[0]>freArr1[0]&&
					freArr2[0]-freArr1[0]<10&&
					freArr2[1]>freArr1[1]&&
					freArr2[1]-freArr1[1]<25){

					flag=1;
				}
			}
		}

		if(!flag&&
			length1>2&&length2>2){

			// fre asc
			__vcorrsort1(freArr1,dbArr1,NULL,NULL,3,0);
			__vcorrsort1(freArr2,dbArr2,NULL,NULL,3,0);

			flag=1;
			for(int i=0;i<2;i++){
				int k=0;

				k=util_calRangeTimes(freArr1[i], freArr2[i], NULL);
				if(k!=1){
					flag=0;
					break;
				}
			}
		}
	}
	else{ // valid ->196+7
		if(length1>10&&length2>10){
			// dB desc
			__vcorrsort1(dbArr1,freArr1,NULL,NULL,length1,1);
			__vcorrsort1(dbArr2,freArr2,NULL,NULL,length2,1);

			if(freArr1[0]>190&&freArr1[0]<204&&
				freArr2[0]>190&&freArr2[0]<204){

				// db1=dbArr1[0];
				// db2=dbArr2[0];
				// if(fabsf(db1-db2)>5.6){
				// 	return 0;
				// }

				// fre asc
				__vcorrsort1(freArr1,dbArr1,NULL,NULL,2,0);
				__vcorrsort1(freArr2,dbArr2,NULL,NULL,2,0);

				flag=1;
				for(int i=0;i<2;i++){
					int k=0;

					k=util_calRangeTimes(freArr1[i], freArr2[i], NULL);
					if(k!=1){
						flag=0;
						break;
					}
				}
			}
		}
		else{ // 330 ->is equal
			// if(length1>4&&length2>4){
			// 	flag=1;
			// 	for(int i=0;i<4;i++){
			// 		int k=0;

			// 		k=util_calRangeTimes(freArr1[i], freArr2[i], NULL);
			// 		if(k!=1){
			// 			flag=0;
			// 			break;
			// 		}
			// 	}
			// }
		}
	}

	return flag;
}

// p 1 isPostive 1 isExp 0 type 0 sum 1 mean
static float __calFlux(float *curArr,float *preArr,int length,int p,int isPostive,int isExp,int type){
	float value=0;

	for(int i=0;i<length;i++){
		float v1=0;

		v1=curArr[i]-preArr[i];
		if(isPostive){
			v1=(v1>0?v1:0);
		}
		else{
			v1=fabsf(v1);
		}

		if(p==2.0){
			v1*=v1;
		}
		else{
			v1=powf(v1, p);
		}

		value+=v1;
	}

	if(type){ // mean
		value/=length;
	}

	if(isExp){
		value=powf(value, 1.0/p); 
	}

	return value;
}

static int __arr_less(float *arr1,int length1,float value,float *arr2){
	int length2=0;

	for(int i=0;i<length1;i++){
		if(arr1[i]<value){
			arr2[length2]=arr1[i];
			length2++;
		}
	}

	return length2;
}

static int __arr_maxIndex(float *arr,int length){
	int index=0;

	float value=0;

	if(!length){
		return 0;
	}

	value=arr[0];
	for(int i=1;i<length;i++){
		if(value<arr[i]){
			value=arr[i];
			index=i;
		}
	}

	return index;
}







