// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/flux_correct.h"
#include "../dsp/fft_algorithm.h"

#include "../stft_algorithm.h"

#include "../classic/trist.h"

#include "_queue.h"
#include "_pitch_ffp.h"

struct OpaquePitchFFP{
	int isContinue;

	STFTObj stftObj;

	int fftLength;
	int slideLength;

	int cutLength; // 12
	int peakLength; // (maxIndex-minIndex+1)/2+1

	int minIndex; // min/maxFre
	int maxIndex;

	int timeLength;

	// stft result ->timeLength*(maxIndex-minIndex+1)
	float *mPowerArr;
	float *mDbArr;

	// timeLength*peakLength
	float *mPeakDbArr;
	float *mPeakFreArr;
	float *mPeakHeightArr;
	int *mIndexArr;
	int *lenArr;
	int *lowFlagArr;

	// filter ->timeLength*peakLength 
	float *mFilterDbArr1; // height-filter
	float *mFilterFreArr1; 
	float *mFilterHeightArr1;
	int *mIndexArr1; 
	int *lenArr1;

	float *mFilterDbArr2; // near-filter
	float *mFilterFreArr2; 
	float *mFilterHeightArr2;
	int *mIndexArr2; 
	int *lenArr2;

	float *mFilterDbArr3; // dB-filter
	float *mFilterFreArr3; 
	float *mFilterHeightArr3;
	int *mIndexArr3;
	int *lenArr3;

	float *maxDBArr; // timeLength
	int *sucessTypeArr; // 1 normal, 2 fast
	float *lightArr;

	float *avgTempArr;
	float *maxTempArr;
	float *percentTempArr;

	// format
	int *formatFlagArr;

	float *formatFreArr1;
	float *formatFreArr2;
	float *formatFreArr3;

	float *formatDbArr1;
	float *formatDbArr2;
	float *formatDbArr3;

	// fast ->timeLength*peakLength 
	float *mFastDbArr2;
	float *mFastFreArr2;
	float *mFastHeightArr2;
	int *mFastIndexArr2;
	int *fastLenArr2;

	float *mFastDbArr3; // dB-filter
	float *mFastFreArr3;
	float *mFastHeightArr3;
	int *mFastIndexArr3;
	int *fastLenArr3;

	float *mFastDbArr4; // cut-filter
	float *mFastFreArr4;
	float *mFastHeightArr4;
	int *mFastIndexArr4;
	int *fastLenArr4;

	// cache
	float *mRealArr; // timeLength*fftLength
	float *mImageArr;

	int *domIndexArr;
	int domLength;

	int samplate;
	WindowType winType;
	float tempBase; // -18 >-36&&<0

	int isDebug;
};

static void __pitchFFPObj_dealData(PitchFFPObj pitchFFPObj,int dataLength);

static void __pitchFFPObj_stft(PitchFFPObj pitchFFPObj,float *dataArr,int dataLength,float *dbArr);
static void __pitchFFPObj_sub(PitchFFPObj pitchFFPObj,float *freArr);

static void __pitchFFPObj_filter(PitchFFPObj pitchFFPObj);
static void __pitchFFPObj_filterFast(PitchFFPObj pitchFFPObj);

static void __pitchFFPObj_filterHeight(PitchFFPObj pitchFFPObj);
static void __pitchFFPObj_filterNear(PitchFFPObj pitchFFPObj);
static void __pitchFFPObj_filterDB(PitchFFPObj pitchFFPObj);
static void __pitchFFPObj_filterRelation(PitchFFPObj pitchFFPObj);

static void __pitchFFPObj_filterFastDB(PitchFFPObj pitchFFPObj);
static void __pitchFFPObj_filterFastCut(PitchFFPObj pitchFFPObj);

static int __pitchFFPObj_preprocess(PitchFFPObj pitchFFPObj,int index);

static void __pitchFFPObj_temporal(PitchFFPObj pitchFFPObj,float *dataArr,int dataLength);

static int __isLowFre(float *freArr,float *heightArr,int *indexArr,int length);
static float __isLight(float *dataArr,int length);
static float __temproal(float *dataArr,int length,float base,float *avgValue,float *percentValue);

static int __arr_rectify(float *dbArr,float *freArr,float *heightArr,int *indexArr,int length);

static int __arr_maxIndex(float *arr,int length);
static int __arr_has(int *arr,int length,int valud);

/***
	samplate 32000
	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int pitchFFPObj_new(PitchFFPObj *pitchFFPObj,
					int *samplate,float *lowFre,float *highFre,
					int *radix2Exp,int *slideLength,WindowType *windowType,
					int *isContinue){
	int status=0;

	int _samplate=32000;
	float _lowFre=27; // 27.5
	float _highFre=4000; // 2093/4186
	int _radix2Exp=12;
	int _slideLength=0;
	WindowType _winType=Window_Hamm;
	int _isContinue=0;

	int fftLength=0;
	int cutLength=12;
	int peakLength=0;

	int minIndex=0; // min/maxFre
	int maxIndex=0;

	STFTObj stftObj=NULL;
	PitchFFPObj pitch=NULL;

	pitch=*pitchFFPObj=(PitchFFPObj )calloc(1,sizeof(struct OpaquePitchFFP ));

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
			_highFre=4000;
		}
	}

	if(radix2Exp){
		if(*radix2Exp>=1&&*radix2Exp<=30){
			_radix2Exp=*radix2Exp;
		}
	}

	if(windowType){
		if(*windowType<=Window_Hamm){
			_winType=*windowType;
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

	minIndex=floorf(_lowFre*fftLength/_samplate);
	maxIndex=ceilf(_highFre*fftLength/_samplate);
	if(maxIndex>=fftLength/2){
		maxIndex=fftLength/2-1;
	}

	if(minIndex>=maxIndex){
		minIndex=3;
		maxIndex=ceilf(4000*fftLength/_samplate);
	}

	peakLength=(maxIndex-minIndex)/2+1;

	stftObj_new(&stftObj, _radix2Exp, &_winType, &_slideLength, &_isContinue);
	
	pitch->isContinue=_isContinue;

	pitch->stftObj=stftObj;

	pitch->fftLength=fftLength;
	pitch->slideLength=_slideLength;

	pitch->cutLength=cutLength;
	pitch->peakLength=peakLength;

	pitch->minIndex=minIndex;
	pitch->maxIndex=maxIndex;

	pitch->samplate=_samplate;
	pitch->winType=_winType;

	pitch->domIndexArr=__vnewi(10, NULL);
	
	return status;
}

void pitchFFPObj_setTempBase(PitchFFPObj pitchFFPObj,float tempBase){

	if(tempBase<0&&tempBase>-36){
		pitchFFPObj->tempBase=tempBase;
	}
}

int pitchFFPObj_calTimeLength(PitchFFPObj pitchFFPObj,int dataLength){
	int timeLen=0;

	timeLen=stftObj_calTimeLength(pitchFFPObj->stftObj, dataLength);
	return timeLen;
}

void pitchFFPObj_pitch(PitchFFPObj pitchFFPObj,float *dataArr,int dataLength,
					float *freArr,float *dbArr){

	// 1. dealData
	__pitchFFPObj_dealData(pitchFFPObj,dataLength);

	// 2. temporal
	__pitchFFPObj_temporal(pitchFFPObj,dataArr,dataLength);

	// 3. stft
	__pitchFFPObj_stft(pitchFFPObj, dataArr, dataLength,dbArr);

	// 4. filter
	__pitchFFPObj_filter(pitchFFPObj);
	__pitchFFPObj_filterFast(pitchFFPObj);

	// 5. sub
	__pitchFFPObj_sub(pitchFFPObj,freArr);

}

int pitchFFPObj_getCorrData(PitchFFPObj pitchFFPObj,float **mCorrArr,float **mDbArr,float **mHeightArr,int **lenArr){
	int mLen=0;

	mLen=pitchFFPObj->peakLength;
	if(mCorrArr){
		*mCorrArr=pitchFFPObj->mFilterFreArr3;
	}

	if(mDbArr){
		*mDbArr=pitchFFPObj->mFilterDbArr3;
	}

	if(mHeightArr){
		*mHeightArr=pitchFFPObj->mFilterHeightArr3;
	}

	if(lenArr){
		*lenArr=pitchFFPObj->lenArr3;
	}

	return mLen;
}

int pitchFFPObj_getCutData(PitchFFPObj pitchFFPObj,float **mCorrArr,float **mDbArr,float **mHeightArr,int **lenArr){
	int mLen=0;

	mLen=pitchFFPObj->peakLength;
	if(mCorrArr){
		*mCorrArr=pitchFFPObj->mFastFreArr4;
	}

	if(mDbArr){
		*mDbArr=pitchFFPObj->mFastDbArr4;
	}

	if(mHeightArr){
		*mHeightArr=pitchFFPObj->mFastHeightArr4;
	}

	if(lenArr){
		*lenArr=pitchFFPObj->fastLenArr4;
	}

	return mLen;
}

int pitchFFPObj_getFlagData(PitchFFPObj pitchFFPObj,int **flagArr){
	int timeLen=0;

	timeLen=pitchFFPObj->timeLength;
	if(flagArr){
		*flagArr=pitchFFPObj->sucessTypeArr;
	}

	return timeLen;
}

int pitchFFPObj_getLightData(PitchFFPObj pitchFFPObj,float **lightArr){
	int timeLen=0;

	timeLen=pitchFFPObj->timeLength;
	if(lightArr){
		*lightArr=pitchFFPObj->lightArr;
	}

	return timeLen;
}

int pitchFFPObj_getFormatData(PitchFFPObj pitchFFPObj,
							int **formatFlagArr,
							float **formatFreArr1,float **formatFreArr2,float **formatFreArr3,
							float **formatDbArr1,float **formatDbArr2,float **formatDbArr3){
	int timeLen=0;

	if(formatFlagArr){
		*formatFlagArr=pitchFFPObj->formatFlagArr;
	}

	if(formatFreArr1){
		*formatFreArr1=pitchFFPObj->formatFreArr1;
	}
	if(formatFreArr2){
		*formatFreArr2=pitchFFPObj->formatFreArr2;
	}
	if(formatFreArr3){
		*formatFreArr3=pitchFFPObj->formatFreArr3;
	}

	if(formatDbArr1){
		*formatDbArr1=pitchFFPObj->formatDbArr1;
	}
	if(formatDbArr2){
		*formatDbArr2=pitchFFPObj->formatDbArr2;
	}
	if(formatDbArr3){
		*formatDbArr3=pitchFFPObj->formatDbArr3;
	}

	return timeLen;
}

int pitchFFPObj_getTemporalData(PitchFFPObj pitchFFPObj,float **avgTempArr,float **maxTempArr,float **percentTempArr){
	int timeLen=0;

	timeLen=pitchFFPObj->timeLength;
	if(avgTempArr){
		*avgTempArr=pitchFFPObj->avgTempArr;
	}

	if(maxTempArr){
		*maxTempArr=pitchFFPObj->maxTempArr;
	}

	if(percentTempArr){
		*percentTempArr=pitchFFPObj->percentTempArr;
	}

	return timeLen;
}

static void __pitchFFPObj_sub(PitchFFPObj pitchFFPObj,float *freArr){
	float *mFilterDbArr3=NULL; // dB-filter
	float *mFilterFreArr3=NULL; 
	float *mFilterHeightArr3=NULL;
	int *lenArr3=NULL;

	float *mFastDbArr3=NULL; // dB-filter
	float *mFastFreArr3=NULL; 
	float *mFastHeightArr3=NULL;
	int *fastLenArr3=NULL;

	float *mFastDbArr4=NULL; // cut-filter
	float *mFastFreArr4=NULL; 
	float *mFastHeightArr4=NULL;
	int *fastLenArr4=NULL;

	int peakLength=0;
	int timeLength=0;
	
	float *corrArr1=NULL;
	float *dbArr1=NULL;
	float *heightArr1=NULL;
	int len1=0;

	float *corrArr2=NULL;
	float *dbArr2=NULL;
	float *heightArr2=NULL;
	int len2=0;

	float *corrArr3=NULL;
	float *dbArr3=NULL;
	float *heightArr3=NULL;
	int len3=0;

	float *lightArr=NULL;

	int flag=0;

	mFilterDbArr3=pitchFFPObj->mFilterDbArr3;
	mFilterFreArr3=pitchFFPObj->mFilterFreArr3;
	mFilterHeightArr3=pitchFFPObj->mFilterHeightArr3;
	lenArr3=pitchFFPObj->lenArr3;

	mFastDbArr3=pitchFFPObj->mFastDbArr3;
	mFastFreArr3=pitchFFPObj->mFastFreArr3;
	mFastHeightArr3=pitchFFPObj->mFastHeightArr3;
	fastLenArr3=pitchFFPObj->fastLenArr3;

	mFastDbArr4=pitchFFPObj->mFastDbArr4;
	mFastFreArr4=pitchFFPObj->mFastFreArr4;
	mFastHeightArr4=pitchFFPObj->mFastHeightArr4;
	fastLenArr4=pitchFFPObj->fastLenArr4;

	peakLength=pitchFFPObj->peakLength;
	timeLength=pitchFFPObj->timeLength;

	lightArr=pitchFFPObj->lightArr;
	for(int i=0;i<timeLength;i++){
		corrArr1=mFilterFreArr3+i*peakLength;
		dbArr1=mFilterDbArr3+i*peakLength;
		heightArr1=mFilterHeightArr3+i*peakLength;
		len1=lenArr3[i];

		corrArr2=mFastFreArr3+i*peakLength;
		dbArr2=mFastDbArr3+i*peakLength;
		heightArr2=mFastHeightArr3+i*peakLength;
		len2=fastLenArr3[i];

		corrArr3=mFastFreArr4+i*peakLength;
		dbArr3=mFastDbArr4+i*peakLength;
		heightArr3=mFastHeightArr4+i*peakLength;
		len3=fastLenArr4[i];

		pitchFFPObj->formatFlagArr[i]=0; // reset

		// load trist!!!

		pitchFFPObj->sucessTypeArr[i]=flag;
	}
}

/***
	1. first 1&2
	2. all-3/4/5
	3. max-1/2/3
****/
static int __pitchFFPObj_preprocess(PitchFFPObj pitchFFPObj,int index){
	float *dbArr=NULL;
	float *freArr=NULL;
	float *hegihtArr=NULL;
	int *indexArr=NULL;

	int len=0;
	int refLen=0;

	int *domIndexArr=NULL;
	int domLen=0;

	float fre1=0;
	float fre2=0;

	float fre3=0;
	float fre4=0;

	float fre5=0;

	int index1=0;
	int index2=0;

	int index3=0;
	int index4=0;

	int index5=0;

	int k1=0;
	int k2=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int vs1=0,vs2=0;
	int vk1=0,vk2=0,vk3=0;

	float _fre=0;
	int _index=0;

	float minHeight=15;

	float *_dbArr=NULL;
	float *_freArr=NULL;
	float *_hegihtArr=NULL;
	int *_indexArr=NULL;

	int _len=0;
	int _offset=0;

	refLen=pitchFFPObj->lenArr3[index];
	domIndexArr=pitchFFPObj->domIndexArr;

	_len=pitchFFPObj->lenArr[index];

	_dbArr=pitchFFPObj->mPeakDbArr+index*pitchFFPObj->peakLength;
	_freArr=pitchFFPObj->mPeakFreArr+index*pitchFFPObj->peakLength;
	_hegihtArr=pitchFFPObj->mPeakHeightArr+index*pitchFFPObj->peakLength;
	_indexArr=pitchFFPObj->mIndexArr+index*pitchFFPObj->peakLength;

	dbArr=__vnew(_len, NULL);
	freArr=__vnew(_len, NULL);
	hegihtArr=__vnew(_len, NULL);
	indexArr=__vnewi(_len, NULL);

	// dB desc
	__vcorrsort1(_dbArr,_freArr,_hegihtArr,_indexArr,_len, 1);

	if(fabsf(_freArr[0]-_freArr[1])>30){
		_offset=0;
		len=_len;
	}
	else{
		_offset=1;
		len=_len-1;
	}

	dbArr[0]=_dbArr[0];
	freArr[0]=_freArr[0];
	hegihtArr[0]=_hegihtArr[0];
	indexArr[0]=_indexArr[0];

	for(int i=1;i+_offset<_len;i++){
		dbArr[i]=_dbArr[i+_offset];
		freArr[i]=_freArr[i+_offset];
		hegihtArr[i]=_hegihtArr[i+_offset];
		indexArr[i]=_indexArr[i+_offset];
	}

	index1=indexArr[0];
	index2=indexArr[1];

	fre1=freArr[0];
	fre2=freArr[1];

	index3=indexArr[2];
	index4=indexArr[3];

	index5=indexArr[4];

	fre3=freArr[2];
	fre4=freArr[3];

	fre5=freArr[4];

	domIndexArr[domLen]=index1;
	domLen++;

	domIndexArr[domLen]=index2;
	domLen++;

	// if(fabsf(fre1-fre2)>30){
	// 	domIndexArr[domLen]=index2;
	// }
	// else{
	// 	domIndexArr[domLen]=index3;
	// }
	// domLen++;

	if(index1>index2){
		_fre=fre1;
		fre1=fre2;
		fre2=_fre;

		_index=index1;
		index1=index2;
		index2=_index;
	}

	__queue_fre2(fre1, fre2, &k1, &k2);
	if(k1==1&&k2==2&&
		(fabsf(fre1*2-fre2)<5||
			(fre1>100&&fre1<120&&
				fabsf(fre1*2-fre2)<15)||
					(fre1>140&&fre1<155&&
						fabsf(fre1*2-fre2)<10))){

		/***
			string-5,1-24-5/7, 2>4>1>5/7, 
					,1234-5
			110+10 ->100~120 ->200~240
		****/
		if(fre3>100&&fre3<120&&
			index3<index1&&index3<index2){

			__queue_fre2(fre3, fre1, &k1, &k2);
			if(k1==1&&k2==2&&
				fabsf(fre3*2-fre1)<4){

				domIndexArr[domLen]=indexArr[2];
				domLen++;

				if(fre4>fre2&&
					hegihtArr[3]>12&&
					(fabsf(fre3*5-fre4)<5||
					fabsf(fre3*7-fre4)<5)){

					domIndexArr[domLen]=indexArr[3];
					domLen++;
				}
			}
		}
		else if(fre1>100&&fre1<120){
			int _count=0;

			for(int i=2;i<len;i++){
				if(freArr[i]>fre2){
					__queue_fre2(fre2/2, freArr[i], &k1, &k2);
					if(k1==1&&(k2==3||k2==4||k2==5)&&
						fabsf(fre2/2*k2-freArr[i])<5){

						domIndexArr[domLen]=indexArr[i];
						domLen++;
					}

					_count++;
					if(_count>=3){
						break;
					}
				}
			}
		}

		// ->236, low65~75, 3>6>2>N, x236/236x
		if(index3<index1&&indexArr[0]<indexArr[1]&&
			((hegihtArr[0]>minHeight&&
				hegihtArr[1]>minHeight)||
			(hegihtArr[0]>minHeight+3&&
				hegihtArr[1]>minHeight-2))){

			if(fre3<150&&fre3>130){
				__queue_fre2(fre3, fre1, &k1, &k2);
				if(k1==2&&k2==3&&
					fabsf(fre3/k1*k2-fre1)<5){

					if(refLen>=3){
						domIndexArr[domLen]=indexArr[2];
						domLen++;
					}
				}
			}
			else{
				if(index4<index1&&indexArr[0]<indexArr[1]&&
					index4>index3&&
					dbArr[2]-dbArr[3]<2){

					if(fre3>90&&fre3<110&&
						fre4<150&&fre4>130){

						__queue_fre2(fre4, fre1, &k1, &k2);
						if(k1==2&&k2==3&&
							fabsf(fre4/k1*k2-fre1)<5){

							if(refLen>=3){
								domIndexArr[domLen]=indexArr[3];
								domLen++;
							}
						}
					}
				}
			}
		}

		// ->234, low50~60, 2>4>3>N, x234/234x, 2>13.1,34>minFre
		if(fre1>100&&fre1<120&&
			((hegihtArr[0]>minHeight&&
				hegihtArr[1]>minHeight)||
			(hegihtArr[0]>minHeight+3&&
				hegihtArr[1]>minHeight-2))){

			int vFlag=1;
			int _count=0;

			// if(freArr[0]>108&&
			// 	freArr[0]<112&&
			// 	dbArr[0]-dbArr[1]>18){

			// 	vFlag=0;
			// }

			if(vFlag){
				for(int i=2;i<len;i++){
					if(freArr[i]>150&&freArr[i]<180&&
						indexArr[i]>index1&&indexArr[i]<index2){

						__queue_fre2(fre1, freArr[i], &k1, &k2);
						if(k1==2&&k2==3&&
							fabsf(fre1/k1*k2-freArr[i])<5){

							if(refLen>=3){
								domIndexArr[domLen]=indexArr[i];
								domLen++;
							}
						}

						_count++;
						if(_count>=3){
							break;
						}
					}
				}
			}
		}

		// ->123, 147+7 
		if(fre1>140&&fre1<154&&
			indexArr[0]>indexArr[1]){

			int _count=0;

			for(int i=2;i<len;i++){
				if(freArr[i]>fre2){
					__queue_fre2(fre1, freArr[i], &k1, &k2);
					if(k1==1&&(k2==3||k2==4)&&
						(fabsf(fre1*k2-freArr[i])<5||
							fabsf(fre1-freArr[i]/k2)<3)){

						domIndexArr[domLen]=indexArr[i];
						domLen++;
					}

					_count++;
					if(_count>=3){
						break;
					}
				}
			}
		}

		// ->234, 80 ->75~90,
		if(fre1>75&&fre1<90){
			for(int i=2;i<len;i++){
				if(freArr[i]>fre2){
					__queue_fre2(fre2, freArr[i], &k1, &k2);
					if(k1==2&&k2==3){
						domIndexArr[domLen]=indexArr[i];
						domLen++;
					}

					break;
				}
			}
		}

		// ->234, low50~60
		// if(fre1>100&&fre1<120){
		// 	for(int i=2;i<len;i++){
		// 		if(freArr[i]<fre1){
		// 			__queue_fre2(freArr[i],fre1,&k1, &k2);
		// 			if(k1==1&&k2==2&&
		// 				hegihtArr[i]>minHeight&&
		// 				fabsf(freArr[i]/k1*k2-fre1)<5){

		// 				if(refLen>=3){
		// 					domIndexArr[domLen]=indexArr[i];
		// 					domLen++;
		// 				}
		// 			}

		// 			break;
		// 		}
		// 	}
		// }
	}
	
	__queue_fre2(fre1, fre2, &k1, &k2);
	if(k1==1&&k2==3&&
		fabsf(fre1*3-fre2)<5){

		/***
			string-6,1267, 2>6>1>7, 2-6>15
			80 ->70~85 ->140~170
			147+5->152, format>85
		****/
		if(freArr[0]>140&&freArr[0]<170&&
			indexArr[0]<indexArr[1]){

			__queue_fre2(fre3, freArr[0], &k1, &k2);
			if(k1==1&&k2==2&&
				fabsf(fre3*2-freArr[0])<4){

				if(hegihtArr[0]>minHeight&&hegihtArr[1]>minHeight){
					domIndexArr[domLen]=indexArr[2];
					domLen++;
				}
			}
		}	
	}

	/***
		string-6,x23x
		80+5 ->75~85
	****/
	if(freArr[0]>150&&freArr[0]<170&&
		indexArr[0]>indexArr[1]){

		int _count=0;

		for(int i=2;i<len;i++){
			if(freArr[i]>freArr[0]){
				__queue_fre2(freArr[0]/2, freArr[i], &k1, &k2);
				if(k1==1&&k2==3&&
					((fabsf(freArr[0]/2*k2-freArr[i])<4)||
					(i==2&&fabsf(freArr[0]/2*k2-freArr[i])<5))){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=3){
					break;
				}
			}
		}
	}

	__queue_fre2(fre1, fre2, &k1, &k2);
	if(freArr[0]>150&&freArr[0]<170&&
		k1==2&&k2==3&&
		fre3<freArr[0]){

		__queue_fre2(fre3, freArr[0], &k1, &k2);
		if(k1==1&&k2==2&&
			fabsf(fre3*2-freArr[0])<4){

			domIndexArr[domLen]=indexArr[2];
			domLen++;
		}
		else if(fre4<freArr[0]&&
			dbArr[2]-dbArr[3]<3&&
			indexArr[2]-indexArr[3]<=3){

			__queue_fre2(fre4, freArr[0], &k1, &k2);
			if(fabsf(fre4*2-freArr[0])<4){
				domIndexArr[domLen]=indexArr[3];
				domLen++;
			}
		}
	}

	__queue_fre2(fre1, fre2, &k1, &k2);
	if(freArr[0]>150&&freArr[0]<170&&
		k1==1&&k2==3){

		if(fre3>freArr[0]){
			__queue_fre2(freArr[0], fre3, &k1, &k2);
			if(k1==2&&k2==3&&
				fabsf(freArr[0]/2-fre3/3)<5){ // error<=15

				domIndexArr[domLen]=indexArr[2];
				domLen++;
			}
		}
		else{
			__queue_fre2(fre3, freArr[0], &k1, &k2);
			if(k1==1&&k2==2&&
				fabsf(fre3*2-freArr[0])<5){

				domIndexArr[domLen]=indexArr[2];
				domLen++;
			}
		}
	}

	// __queue_fre2(fre1, fre2, &k1, &k2);
	if(freArr[0]>150&&freArr[0]<170){ // &&k1==1&&k2==2

		int _count=0;

		for(int i=1;i<len;i++){
			if(freArr[i]>freArr[0]){
				// __queue_fre2(freArr[0], freArr[i], &k1, &k2);
				if(fabsf(freArr[0]/2-freArr[i]/3)<5){ // error<=15 k1==2&&k2==3&&

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=3){
					break;
				}
			}
		}
	}

	/***
		string-5,1x23
		110+10 ->100~120, 1234/134,3>1>2>4
	****/
	if(freArr[0]>100&&freArr[0]<120&&
		indexArr[0]<indexArr[1]&&
		refLen>3){

		int _count=0;

		for(int i=1;i<len;i++){
			if(freArr[i]>freArr[0]){
				__queue_fre2(freArr[0], freArr[i], &k1, &k2);
				if(k1==1&&(k2==2||k2==3||k2==4)){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=3){
					break;
				}
			}
		}
	}

	__queue_fre2(fre1, fre2, &k1, &k2);
	if(freArr[0]/2>100&&freArr[0]/2<120&&
		indexArr[0]>indexArr[1]&&
		k1==1&&k2==2){

		int _count=0;

		for(int i=1;i<len;i++){
			if(freArr[i]>freArr[0]){
				__queue_fre2(freArr[0]/2, freArr[i], &k1, &k2);
				if(k1==1&&(k2==3||k2==4)){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=2){
					break;
				}
			}
		}
	}

	__queue_fre2(fre1, fre2, &k1, &k2);
	if(freArr[0]>315&&freArr[0]<345&&
		indexArr[0]>indexArr[1]&&
		k1==1&&k2==3){ // dbArr[0]-dbArr[1]<4||pitchFFPObj->lightArr[index]>0.98

		int _count=0;

		for(int i=1;i<len;i++){
			if(freArr[i]>freArr[0]){
				__queue_fre2(fre1, freArr[i], &k1, &k2);
				if(k1==1&&k2==4){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=1){
					break;
				}
			}
		}
	}

	if(freArr[0]>200&&freArr[0]<240&&
		indexArr[0]<indexArr[1]&&
		k1==1&&k2==2&&
		fabsf(fre1*2-fre2)<5){ // 24, ->1234, 2>4

		int _count=0;

		for(int i=2;i<len;i++){
			if(freArr[i]>freArr[0]){
				__queue_fre2(fre1, freArr[i], &k1, &k2);
				if(k1==2&&k2==3){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=2){
					break;
				}
			}
		}
	}

	if(freArr[0]>200&&freArr[0]<240&&
		indexArr[0]>indexArr[1]&&
		k1==1&&k2==2&&
		fabsf(fre1*2-fre2)<5){ // 12, ->1234, 2>1

		int _count=0;

		for(int i=2;i<len;i++){
			if(freArr[i]>freArr[0]){
				__queue_fre2(fre1, freArr[i], &k1, &k2);
				if(k1==2&&k2==3){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=2){
					break;
				}
			}
		}
	}

	/***
		string-4, x123/x136/x1x2, 1-2/3>18
		147+7 ->140~154
	****/
	if(freArr[0]>140&&freArr[0]<154&&
		indexArr[0]>indexArr[1]){

		int _count=0;

		for(int i=2;i<len;i++){
			if(freArr[i]>freArr[0]){
				__queue_fre2(freArr[0], freArr[i], &k1, &k2);
				if(k1==1&&(k2==2||k2==3)&&
					fabsf(freArr[0]*k2-freArr[i])<5){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=3){
					break;
				}
			}
		}
	}

	if(freArr[0]>280&&freArr[0]<310){
		int _count=0;

		for(int i=1;i<len;i++){
			if(freArr[i]<freArr[0]){
				__queue_fre2(freArr[i],freArr[0], &k1, &k2);
				if(k1==1&&k2==2&&
					fabsf(freArr[i]*k2-freArr[0])<8){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=2){
					break;
				}
			}
		}

		_count=0;
		for(int i=1;i<len;i++){
			if(freArr[i]>freArr[0]){
				__queue_fre2(freArr[0]/2, freArr[i],&k1, &k2);
				if(k1==1&&(k2==3||k2==4)&&
					fabsf(freArr[0]/2*k2-freArr[i])<5){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=3){
					break;
				}
			}
		}
	}

	/***
		string-3, x13x, 1-3>24 || x-1<3&&1-3>12
		197+7 ->190~204
	****/
	if((freArr[0]>190&&freArr[0]<204&&
		indexArr[0]>indexArr[1])||
		(freArr[1]>190&&freArr[0]<204&&
			dbArr[0]-dbArr[1]<3&&
			indexArr[0]<indexArr[1])){

		int _count=0;
		float _fre=0;

		if(freArr[0]>190&&freArr[0]<204){
			_fre=freArr[0];
		}
		else{
			_fre=freArr[1];
		}

		for(int i=2;i<len;i++){
			if(freArr[i]>_fre){
				__queue_fre2(_fre, freArr[i], &k1, &k2);
				if(k1==1&&(k2==2||k2==3)&&
					fabsf(_fre*k2-freArr[i])<4){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=3){
					break;
				}
			}
		}
	}

	/***
		string-2,123
		247 -> >220
	****/
	__queue_fre2(fre1, fre2, &k1, &k2);
	if(freArr[0]>220&&
		indexArr[0]<indexArr[1]&&
		k1==1&&k2==2&&
		fabsf(fre1*2-fre2)<5&&
		refLen>3){

		int _count=0;

		for(int i=2;i<len;i++){
			if(freArr[i]>freArr[1]){
				__queue_fre2(freArr[0], freArr[i], &k1, &k2);
				if(k1==1&&k2==3){

					domIndexArr[domLen]=indexArr[i];
					domLen++;
				}

				_count++;
				if(_count>=2){
					break;
				}
			}
		}
	}

	// fre asc
	__vcorrsort1(_freArr,_dbArr,_hegihtArr,_indexArr,_len, 0);

	free(dbArr);
	free(freArr);
	free(hegihtArr);
	free(indexArr);

	pitchFFPObj->domLength=domLen;
	return domLen;
}

static void __pitchFFPObj_filterFast(PitchFFPObj pitchFFPObj){
	int timeLength=0;
	int peakLength=0; 

	float *mPeakDbArr=NULL;
	float *mPeakFreArr=NULL;
	float *mPeakHeightArr=NULL;
	int *mIndexArr=NULL;
	int *lenArr=NULL;

	float *mFastDbArr2=NULL;
	float *mFastFreArr2=NULL;
	float *mFastHeightArr2=NULL;
	int *mFastIndexArr2=NULL;
	int *fastLenArr2=NULL;

	int *domIndexArr=NULL;
	int domLen=0;

	int len=0;
	int len1=0;

	float curDb=0;
	float nexDb=0;

	float curFre=0;
	float nexFre=0;

	float minHeight=15;
	float minFre=30;

	timeLength=pitchFFPObj->timeLength;
	peakLength=pitchFFPObj->peakLength;

	mPeakDbArr=pitchFFPObj->mPeakDbArr;
	mPeakFreArr=pitchFFPObj->mPeakFreArr;
	mPeakHeightArr=pitchFFPObj->mPeakHeightArr;
	mIndexArr=pitchFFPObj->mIndexArr;
	lenArr=pitchFFPObj->lenArr;

	mFastDbArr2=pitchFFPObj->mFastDbArr2;
	mFastFreArr2=pitchFFPObj->mFastFreArr2;
	mFastHeightArr2=pitchFFPObj->mFastHeightArr2;
	mFastIndexArr2=pitchFFPObj->mFastIndexArr2;
	fastLenArr2=pitchFFPObj->fastLenArr2;

	domIndexArr=pitchFFPObj->domIndexArr;

	for(int i=0;i<timeLength;i++){
		len1=0;
		len=lenArr[i];

		domLen=__pitchFFPObj_preprocess(pitchFFPObj,i);

		for(int j=0;j<len;j++){
			if(mPeakHeightArr[i*peakLength+j]>minHeight||
				__arr_has(domIndexArr, domLen, mIndexArr[i*peakLength+j])){

				int _index=0;

				curFre=mPeakFreArr[i*peakLength+j];
				curDb=mPeakDbArr[i*peakLength+j];

				nexFre=0;
				nexDb=0;
				for(int k=j+1;k<len;k++){
					if(mPeakHeightArr[i*peakLength+k]>minHeight||
						__arr_has(domIndexArr, domLen, mIndexArr[i*peakLength+k])){

						nexFre=mPeakFreArr[i*peakLength+k];
						nexDb=mPeakDbArr[i*peakLength+k];

						_index=k;
						break;
					}
				}

				if(nexFre){
					if(nexFre-curFre<minFre){ // minFre sucess
						if(curDb<nexDb){
							mFastDbArr2[i*peakLength+len1]=mPeakDbArr[i*peakLength+_index];
							mFastFreArr2[i*peakLength+len1]=mPeakFreArr[i*peakLength+_index];
							mFastHeightArr2[i*peakLength+len1]=mPeakHeightArr[i*peakLength+_index];
							mFastIndexArr2[i*peakLength+len1]=mIndexArr[i*peakLength+_index];

							len1++;
						}
						else{
							mFastDbArr2[i*peakLength+len1]=mPeakDbArr[i*peakLength+j];
							mFastFreArr2[i*peakLength+len1]=mPeakFreArr[i*peakLength+j];
							mFastHeightArr2[i*peakLength+len1]=mPeakHeightArr[i*peakLength+j];
							mFastIndexArr2[i*peakLength+len1]=mIndexArr[i*peakLength+j];

							len1++;
						}

						j=_index; // update
					}
					else{
						mFastDbArr2[i*peakLength+len1]=mPeakDbArr[i*peakLength+j];
						mFastFreArr2[i*peakLength+len1]=mPeakFreArr[i*peakLength+j];
						mFastHeightArr2[i*peakLength+len1]=mPeakHeightArr[i*peakLength+j];
						mFastIndexArr2[i*peakLength+len1]=mIndexArr[i*peakLength+j];

						len1++;

						// mFastDbArr2[i*peakLength+len1]=mPeakDbArr[i*peakLength+_index];
						// mFastFreArr2[i*peakLength+len1]=mPeakFreArr[i*peakLength+_index];
						// mFastHeightArr2[i*peakLength+len1]=mPeakHeightArr[i*peakLength+_index];
						// mFastIndexArr2[i*peakLength+len1]=mIndexArr[i*peakLength+_index];

						// len1++;
					}
				}
				else{ // end
					mFastDbArr2[i*peakLength+len1]=mPeakDbArr[i*peakLength+j];
					mFastFreArr2[i*peakLength+len1]=mPeakFreArr[i*peakLength+j];
					mFastHeightArr2[i*peakLength+len1]=mPeakHeightArr[i*peakLength+j];
					mFastIndexArr2[i*peakLength+len1]=mIndexArr[i*peakLength+j];

					len1++;
				}
			}
		}

		fastLenArr2[i]=len1;
	}

	__pitchFFPObj_filterFastDB(pitchFFPObj);
	__pitchFFPObj_filterFastCut(pitchFFPObj);
}

static void __pitchFFPObj_filter(PitchFFPObj pitchFFPObj){

	__pitchFFPObj_filterHeight(pitchFFPObj);
	__pitchFFPObj_filterNear(pitchFFPObj);
	__pitchFFPObj_filterDB(pitchFFPObj);
	__pitchFFPObj_filterRelation(pitchFFPObj);

}

// minHeight
static void __pitchFFPObj_filterHeight(PitchFFPObj pitchFFPObj){
	int timeLength=0;
	int peakLength=0; 

	float *mPeakDbArr=NULL;
	float *mPeakFreArr=NULL;
	float *mPeakHeightArr=NULL;
	int *mIndexArr=NULL;
	int *lenArr=NULL;
	int *lowFlagArr=NULL;

	float *mFilterDbArr1=NULL; // height-filter
	float *mFilterFreArr1=NULL; 
	float *mFilterHeightArr1=NULL;
	int *mIndexArr1=NULL; 
	int *lenArr1=NULL;

	int len=0;
	int len1=0;

	float minHeight=15;

	timeLength=pitchFFPObj->timeLength;
	peakLength=pitchFFPObj->peakLength;

	mPeakDbArr=pitchFFPObj->mPeakDbArr;
	mPeakFreArr=pitchFFPObj->mPeakFreArr;
	mPeakHeightArr=pitchFFPObj->mPeakHeightArr;
	mIndexArr=pitchFFPObj->mIndexArr;
	lenArr=pitchFFPObj->lenArr;
	lowFlagArr=pitchFFPObj->lowFlagArr;

	mFilterDbArr1=pitchFFPObj->mFilterDbArr1;
	mFilterFreArr1=pitchFFPObj->mFilterFreArr1;
	mFilterHeightArr1=pitchFFPObj->mFilterHeightArr1;
	mIndexArr1=pitchFFPObj->mIndexArr1;
	lenArr1=pitchFFPObj->lenArr1;

	for(int i=0;i<timeLength;i++){
		int start=0;

		int firstIndex=0;
		int secondIndex=0;

		len1=0;
		len=lenArr[i];

		if(len>=2){
			start=2;
			len1=2;
		}
		else if(len>=1){
			start=1;
			len1=1;
		}

		for(int j=0;j<len1;j++){
			mFilterDbArr1[i*peakLength+j]=mPeakDbArr[i*peakLength+j];
			mFilterFreArr1[i*peakLength+j]=mPeakFreArr[i*peakLength+j];
			mFilterHeightArr1[i*peakLength+j]=mPeakHeightArr[i*peakLength+j];
			mIndexArr1[i*peakLength+j]=mIndexArr[i*peakLength+j];

			if(j==0){
				firstIndex=mIndexArr[i*peakLength+j];
			}
			else if(j==1){
				secondIndex=mIndexArr[i*peakLength+j];
			}
		}

		if(lowFlagArr[i]){
			for(int j=start;j<len;j++){
				if(mPeakHeightArr[i*peakLength+j]>minHeight){
					mFilterDbArr1[i*peakLength+len1]=mPeakDbArr[i*peakLength+j];
					mFilterFreArr1[i*peakLength+len1]=mPeakFreArr[i*peakLength+j];
					mFilterHeightArr1[i*peakLength+len1]=mPeakHeightArr[i*peakLength+j];
					mIndexArr1[i*peakLength+len1]=mIndexArr[i*peakLength+j];

					len1++;
				}
			}
		}
		else{
			// fre asc ->start~len1
			__vcorrsort1(mPeakFreArr+(i*peakLength+start),
						mPeakDbArr+(i*peakLength+start), 
						mPeakHeightArr+(i*peakLength+start),
						mIndexArr+(i*peakLength+start),
						len-start, 0);

			for(int j=start;j<len-1;j++){
				if(mPeakHeightArr[i*peakLength+j]>minHeight){
					float curDb=0;
					float preDb=0;
					float nexDb=0;

					float curFre=0;
					float preFre=0;
					float nexFre=0;

					float curHeight=0;
					float preHeight=0;
					float nexHeight=0;

					int curIndex=0;
					int preIndex=0;
					int nexIndex=0;

					int flag=0;

					curDb=mPeakDbArr[i*peakLength+j];
					preDb=mPeakDbArr[i*peakLength+j-1];
					nexDb=mPeakDbArr[i*peakLength+j+1];

					curFre=mPeakFreArr[i*peakLength+j];
					preFre=mPeakFreArr[i*peakLength+j-1];
					nexFre=mPeakFreArr[i*peakLength+j+1];

					curHeight=mPeakHeightArr[i*peakLength+j];
					preHeight=mPeakHeightArr[i*peakLength+j-1];
					nexHeight=mPeakHeightArr[i*peakLength+j+1];

					curIndex=mIndexArr[i*peakLength+j];
					preIndex=mIndexArr[i*peakLength+j-1];
					nexIndex=mIndexArr[i*peakLength+j+1];

					if(firstIndex){
						if(firstIndex>preIndex&&firstIndex<curIndex){
							preHeight=minHeight+1;
						}
					}

					if(secondIndex){
						if(secondIndex>preIndex&&secondIndex<curIndex){
							preHeight=minHeight+1;
						}
					}

					if(firstIndex){
						if(firstIndex>curIndex&&firstIndex<nexIndex){
							nexHeight=minHeight+1;
						}
					}

					if(secondIndex){
						if(secondIndex>curIndex&&secondIndex<nexIndex){
							nexHeight=minHeight+1;
						}
					}

					// if(j==len-1){ // end
					// 	nexHeight=minHeight+1;
					// }

					// if(curIndex-preIndex>=8){
					// 	preHeight=minHeight+1;
					// }

					// if(nexIndex-curIndex>=8){
					// 	nexHeight=minHeight+1;
					// }

					if(curDb>-60){
						if(((curDb-preDb>12)||preHeight>minHeight)&&
							((curDb-nexDb>12)||nexHeight>minHeight)){

							flag=1;
						}
					}
					else{
						float _base=0;

						_base=12;
						if(curHeight>minHeight+4){
							_base=11;
						}

						if(((curDb-preDb>_base)||(preHeight>minHeight&&curIndex-preIndex>3))&&
							((curDb-nexDb>_base)||(nexHeight>minHeight&&nexIndex-curIndex>3))){

							flag=1;
						}
					}

					if(flag){

						mFilterDbArr1[i*peakLength+len1]=mPeakDbArr[i*peakLength+j];
						mFilterFreArr1[i*peakLength+len1]=mPeakFreArr[i*peakLength+j];
						mFilterHeightArr1[i*peakLength+len1]=mPeakHeightArr[i*peakLength+j];
						mIndexArr1[i*peakLength+len1]=mIndexArr[i*peakLength+j];

						len1++;
					}
				}
			}
		}
		
		lenArr1[i]=len1;

		// fre asc 
		__vcorrsort1(mPeakFreArr+i*peakLength,
					mPeakDbArr+i*peakLength, 
					mPeakHeightArr+i*peakLength,
					mIndexArr+i*peakLength,
					len, 0);

		// fre asc
		__vcorrsort1(mFilterFreArr1+i*peakLength,
					mFilterDbArr1+i*peakLength, 
					mFilterHeightArr1+i*peakLength,
					mIndexArr1+i*peakLength,
					len1, 0);
	}
}

// minFre
static void __pitchFFPObj_filterNear(PitchFFPObj pitchFFPObj){
	int timeLength=0;
	int peakLength=0; 

	float *mFilterDbArr1=NULL; // height-filter
	float *mFilterFreArr1=NULL; 
	float *mFilterHeightArr1=NULL;
	int *mIndexArr1=NULL; 
	int *lenArr1=NULL;

	float *mFilterDbArr2=NULL; // near-filter
	float *mFilterFreArr2=NULL; 
	float *mFilterHeightArr2=NULL;
	int *mIndexArr2=NULL; 
	int *lenArr2=NULL;

	int len1=0;
	int len2=0;

	float curDb=0;
	float nexDb=0;
	float nnDb=0;

	float curFre=0;
	float nexFre=0;
	float nnFre=0;

	float curHeight=0;
	float nexHeight=0;

	int lastFlag=1;

	float minFre=30;

	timeLength=pitchFFPObj->timeLength;
	peakLength=pitchFFPObj->peakLength;

	mFilterDbArr1=pitchFFPObj->mFilterDbArr1;
	mFilterFreArr1=pitchFFPObj->mFilterFreArr1;
	mFilterHeightArr1=pitchFFPObj->mFilterHeightArr1;
	mIndexArr1=pitchFFPObj->mIndexArr1;
	lenArr1=pitchFFPObj->lenArr1;

	mFilterDbArr2=pitchFFPObj->mFilterDbArr2;
	mFilterFreArr2=pitchFFPObj->mFilterFreArr2;
	mFilterHeightArr2=pitchFFPObj->mFilterHeightArr2;
	mIndexArr2=pitchFFPObj->mIndexArr2;
	lenArr2=pitchFFPObj->lenArr2;

	for(int i=0;i<timeLength;i++){
		lastFlag=1;
		len2=0;
		len1=lenArr1[i];
		
		for(int j=0;j<len1-1;j++){
			int _index=0;

			curFre=mFilterFreArr1[i*peakLength+j];
			nexFre=mFilterFreArr1[i*peakLength+j+1];

			_index=j;
			if(nexFre-curFre<minFre){
				curDb=mFilterDbArr1[i*peakLength+j];
				nexDb=mFilterDbArr1[i*peakLength+j+1];

				if(j==len1-2){
					lastFlag=0;
				}

				if(curDb<nexDb){
					_index=j+1;
					if(j+2<len1){
						nnFre=mFilterFreArr1[i*peakLength+j+2];
						nnDb=mFilterDbArr1[i*peakLength+j+2];

						if(nnFre-nexFre<minFre&&nexDb>nnDb){
							j++;
						}
					}
				}

				j++;
			}
			else if(nexFre-curFre<2*minFre){ // sub>15||<-70
				curDb=mFilterDbArr1[i*peakLength+j];
				nexDb=mFilterDbArr1[i*peakLength+j+1];

				// if(curDb>nexDb){
				// 	if(curDb-nexDb>15||nexDb<-70){ // sucess
				// 		if(j==len1-2){
				// 			lastFlag=0;
				// 		}

				// 		j++;
				// 	}
				// }
				// else{
				// 	if(nexDb-curDb>15||curDb<-70){ // sucess
				// 		if(j==len1-2){
				// 			lastFlag=0;
				// 		}

				// 		_index=j+1;

				// 		j++;
				// 	}
				// }

				// if(!pitchFFPObj->lowFlagArr[i]){
				// 	if(j==len1-2){
				// 		lastFlag=0;
				// 	}

				// 	if(curDb<nexDb){
				// 		_index=j+1;
				// 		if(j+2<len1){
				// 			nnFre=mFilterFreArr1[i*peakLength+j+2];
				// 			nnDb=mFilterDbArr1[i*peakLength+j+2];

				// 			if(nnFre-nexFre<2*minFre&&nexDb>nnDb){
				// 				j++;
				// 			}
				// 		}
				// 	}

				// 	j++;
				// }
				// else{
				// 	if(curDb>nexDb){
				// 		if(curDb-nexDb>15||nexDb<-70){ // sucess
				// 			if(j==len1-2){
				// 				lastFlag=0;
				// 			}

				// 			j++;
				// 		}
				// 	}
				// 	else{
				// 		if(nexDb-curDb>15||curDb<-70){ // sucess
				// 			if(j==len1-2){
				// 				lastFlag=0;
				// 			}

				// 			_index=j+1;

				// 			j++;
				// 		}
				// 	}
				// }
			}

			mFilterDbArr2[i*peakLength+len2]=mFilterDbArr1[i*peakLength+_index];
			mFilterFreArr2[i*peakLength+len2]=mFilterFreArr1[i*peakLength+_index];
			mFilterHeightArr2[i*peakLength+len2]=mFilterHeightArr1[i*peakLength+_index];
			mIndexArr2[i*peakLength+len2]=mIndexArr1[i*peakLength+_index];

			len2++;
		}

		if(lastFlag){
			mFilterDbArr2[i*peakLength+len2]=mFilterDbArr1[i*peakLength+len1-1];
			mFilterFreArr2[i*peakLength+len2]=mFilterFreArr1[i*peakLength+len1-1];
			mFilterHeightArr2[i*peakLength+len2]=mFilterHeightArr1[i*peakLength+len1-1];
			mIndexArr2[i*peakLength+len2]=mIndexArr1[i*peakLength+len1-1];

			len2++;
		}

		lenArr2[i]=len2;
	}
}

// minDB
static void __pitchFFPObj_filterDB(PitchFFPObj pitchFFPObj){
	int timeLength=0;
	int peakLength=0; 

	float *mFilterDbArr2=NULL; // near-filter
	float *mFilterFreArr2=NULL; 
	float *mFilterHeightArr2=NULL;
	int *mIndexArr2=NULL; 
	int *lenArr2=NULL;

	float *mFilterDbArr3=NULL; // dB-filter
	float *mFilterFreArr3=NULL; 
	float *mFilterHeightArr3=NULL;
	int *mIndexArr3=NULL; 
	int *lenArr3=NULL;

	int len2=0;
	int len3=0;

	float *maxDBArr=NULL;
	float maxDB=0;

	float minDB=14.5; // 15 minDB

	timeLength=pitchFFPObj->timeLength;
	peakLength=pitchFFPObj->peakLength;

	mFilterDbArr2=pitchFFPObj->mFilterDbArr2;
	mFilterFreArr2=pitchFFPObj->mFilterFreArr2;
	mFilterHeightArr2=pitchFFPObj->mFilterHeightArr2;
	mIndexArr2=pitchFFPObj->mIndexArr2;
	lenArr2=pitchFFPObj->lenArr2;

	mFilterDbArr3=pitchFFPObj->mFilterDbArr3;
	mFilterFreArr3=pitchFFPObj->mFilterFreArr3;
	mFilterHeightArr3=pitchFFPObj->mFilterHeightArr3;
	mIndexArr3=pitchFFPObj->mIndexArr3;
	lenArr3=pitchFFPObj->lenArr3;

	maxDBArr=pitchFFPObj->maxDBArr;

	for(int i=0;i<timeLength;i++){
		int start=0;
		int _index=0;

		len3=0;
		len2=lenArr2[i];
		maxDB=maxDBArr[i];

		// -78 filter => len2->lne3
		for(int j=0;j<len2;j++){
			if(mFilterDbArr2[i*peakLength+j]>-100){
				mFilterDbArr3[i*peakLength+len3]=mFilterDbArr2[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len3]=mFilterFreArr2[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len3]=mFilterHeightArr2[i*peakLength+j];
				mIndexArr3[i*peakLength+len3]=mIndexArr2[i*peakLength+j];

				len3++;
			}
		}

		// filter three continue >15+5 =>len3->len2
		{
			float _db1=0;
			float _db2=0;
			float _db3=0;
			float _db4=0;
			float _db5=0;

			len2=0;
			for(int j=0;j<len3;j++){
				mFilterDbArr3[i*peakLength+len2]=mFilterDbArr3[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len2]=mFilterFreArr3[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len2]=mFilterHeightArr3[i*peakLength+j];
				mIndexArr3[i*peakLength+len2]=mIndexArr3[i*peakLength+j];
				len2++;

				if(j+4<len3){
					_db1=mFilterDbArr3[i*peakLength+j];
					_db2=mFilterDbArr3[i*peakLength+j+1];
					_db3=mFilterDbArr3[i*peakLength+j+2];
					_db4=mFilterDbArr3[i*peakLength+j+3];
					_db5=mFilterDbArr3[i*peakLength+j+4];
					if(_db1-_db2>minDB+5&&_db1-_db3>minDB+5&&_db1-_db4>minDB+5&&
						_db5-_db2>minDB+5&&_db5-_db3>minDB+5&&_db5-_db4>minDB+5){ // jump

						j=j+3;
					}
				}
			}
		}

		// filter two continue >15 =>len2->len3
		{
			float _db1=0;
			float _db2=0;
			float _db3=0;
			float _db4=0;

			len3=0;
			for(int j=0;j<len2;j++){
				mFilterDbArr3[i*peakLength+len3]=mFilterDbArr3[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len3]=mFilterFreArr3[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len3]=mFilterHeightArr3[i*peakLength+j];
				mIndexArr3[i*peakLength+len3]=mIndexArr3[i*peakLength+j];
				len3++;

				if(j+3<len2){
					_db1=mFilterDbArr3[i*peakLength+j];
					_db2=mFilterDbArr3[i*peakLength+j+1];
					_db3=mFilterDbArr3[i*peakLength+j+2];
					_db4=mFilterDbArr3[i*peakLength+j+3];
					if(_db1-_db2>minDB&&_db1-_db3>minDB&&
						_db4-_db2>minDB&&_db4-_db3>minDB){ // jump

						j=j+2;
					}
				}
				// else if(j+2<len2){
				// 	_db1=mFilterDbArr3[i*peakLength+j];
				// 	_db2=mFilterDbArr3[i*peakLength+j+1];
				// 	_db3=mFilterDbArr3[i*peakLength+j+2];
				// 	if(_db1-_db2>minDB&&_db1-_db3>minDB){ // last jump
						
				// 		j=j+1;
				// 	}
				// }
			}
		}

		// left -> first <15 => len3->len2
		len2=0;
		_index=__arr_maxIndex(mFilterDbArr3+i*peakLength,len3);
		if(_index>6){ // 240->3520 is large
			_index=0;
		}
		for(int j=0;j<=_index;j++){
			
			if(maxDB-mFilterDbArr3[i*peakLength+j]<minDB||
				mFilterDbArr3[i*peakLength+j]>-42){

				start=j;

				mFilterDbArr3[i*peakLength+len2]=mFilterDbArr3[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len2]=mFilterFreArr3[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len2]=mFilterHeightArr3[i*peakLength+j];
				mIndexArr3[i*peakLength+len2]=mIndexArr3[i*peakLength+j];

				len2++;
			}
		}

		// median -> relative near-filter, not cur dB-filter ???
		for(int j=start+1;j<len3-1;j++){
			if(mFilterDbArr3[i*peakLength+j-1]-mFilterDbArr3[i*peakLength+j]<minDB||
				mFilterDbArr3[i*peakLength+j+1]-mFilterDbArr3[i*peakLength+j]<minDB){

				mFilterDbArr3[i*peakLength+len2]=mFilterDbArr3[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len2]=mFilterFreArr3[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len2]=mFilterHeightArr3[i*peakLength+j];
				mIndexArr3[i*peakLength+len2]=mIndexArr3[i*peakLength+j];

				len2++;
			}
		}

		// end 
		if(len3>1&&start<len3-1){ // for median fre/light
			if(mFilterDbArr3[i*peakLength+len3-2]-mFilterDbArr3[i*peakLength+len3-1]<minDB||
				len3==3||len3==2){

				mFilterDbArr3[i*peakLength+len2]=mFilterDbArr3[i*peakLength+len3-1];
				mFilterFreArr3[i*peakLength+len2]=mFilterFreArr3[i*peakLength+len3-1];
				mFilterHeightArr3[i*peakLength+len2]=mFilterHeightArr3[i*peakLength+len3-1];
				mIndexArr3[i*peakLength+len2]=mIndexArr3[i*peakLength+len3-1];

				len2++;
			}
		}
		
		lenArr3[i]=len2;
	}

	if(pitchFFPObj->isDebug){
		printf("stft start !!!\n");
		for(int i=0;i<timeLength;i++){
			int len=0;

			len=lenArr3[i];

			if(pitchFFPObj->lowFlagArr[i]){
				printf("index[%d], %0.3f: low \n",i,1.0*i*pitchFFPObj->slideLength/pitchFFPObj->samplate);
			}
			else{
				printf("index[%d], %0.3f:\n",i,1.0*i*pitchFFPObj->slideLength/pitchFFPObj->samplate);
			}

			printf("	freArr:\n");
			printf("		");
			for(int j=0;j<len;j++){
				printf("[%d]:%0.1f,",j,mFilterFreArr3[i*peakLength+j]);
			}
			printf("\n");

			printf("	dbArr:\n");
			printf("		");
			for(int j=0;j<len;j++){
				printf("[%d]:%0.1f,",j,mFilterDbArr3[i*peakLength+j]);
			}
			printf("\n");

			printf("	heightArr:\n");
			printf("		");
			for(int j=0;j<len;j++){
				printf("[%d]:%0.1f,",j,mFilterHeightArr3[i*peakLength+j]);
			}
			printf("\n");

			printf("\n");

		}
		printf("stft end !!!\n");
	}
}

static void __pitchFFPObj_filterRelation(PitchFFPObj pitchFFPObj){
	int timeLength=0;
	int peakLength=0; 

	float *mFilterDbArr3=NULL; // dB-filter
	float *mFilterFreArr3=NULL; 
	float *mFilterHeightArr3=NULL;
	int *mIndexArr3=NULL; 
	int *lenArr3=NULL;

	int len3=0;

	int index1=0;
	float minDB=24; // minDB

	timeLength=pitchFFPObj->timeLength;
	peakLength=pitchFFPObj->peakLength;

	mFilterDbArr3=pitchFFPObj->mFilterDbArr3;
	mFilterFreArr3=pitchFFPObj->mFilterFreArr3;
	mFilterHeightArr3=pitchFFPObj->mFilterHeightArr3;
	mIndexArr3=pitchFFPObj->mIndexArr3;
	lenArr3=pitchFFPObj->lenArr3;

	for(int i=0;i<timeLength;i++){
		int flag=0;

		int start=0;
		int end=0;

		len3=lenArr3[i];
		index1=__arr_maxIndex(mFilterDbArr3+i*peakLength, len3);
		if(len3>index1+1&&len3>=12){
			if(index1<=1&&
				mFilterFreArr3[i*peakLength+index1]>190&&
				mFilterFreArr3[i*peakLength+index1]<205){

				int k1=0;
				int k2=0;

				__queue_fre2(mFilterFreArr3[i*peakLength+index1], mFilterFreArr3[i*peakLength+index1+1],
							&k1, &k2);

				if(k1==1&&k2==2){
					start=index1+1;
					for(int j=start+1;j<lenArr3[i];j++){
						if(mFilterDbArr3[i*peakLength+start]-mFilterDbArr3[i*peakLength+j]>minDB){
							__queue_fre2(mFilterFreArr3[i*peakLength+index1], mFilterFreArr3[i*peakLength+j],
									&k1, &k2);

							if(k1==1){
								end=j;
								break;
							}
						}
						else{
							end=j;
							break;
						}
					}
				}

				if(end-start>1&&end-start<4){ // 1||2
					flag=1;
				}
			}
		}

		if(flag){
			for(int j=start+1,k=end;j<len3&&k<len3;j++,k++){
				mFilterDbArr3[i*peakLength+j]=mFilterDbArr3[i*peakLength+k];
				mFilterFreArr3[i*peakLength+j]=mFilterFreArr3[i*peakLength+k];
				mFilterHeightArr3[i*peakLength+j]=mFilterHeightArr3[i*peakLength+k];
				mIndexArr3[i*peakLength+j]=mIndexArr3[i*peakLength+k];
			}

			lenArr3[i]=len3-(end-start)+1;
		}
	}
}

static void __pitchFFPObj_filterFastDB(PitchFFPObj pitchFFPObj){
	int timeLength=0;
	int peakLength=0; 

	float *mFilterDbArr2=NULL; // near-filter
	float *mFilterFreArr2=NULL; 
	float *mFilterHeightArr2=NULL;
	int *mIndexArr2=NULL; 
	int *lenArr2=NULL;

	float *mFilterDbArr3=NULL; // dB-filter
	float *mFilterFreArr3=NULL; 
	float *mFilterHeightArr3=NULL;
	int *mIndexArr3=NULL; 
	int *lenArr3=NULL;

	int len2=0;
	int len3=0;

	float *maxDBArr=NULL;
	float maxDB=0;

	float minDB=15; // minDB

	timeLength=pitchFFPObj->timeLength;
	peakLength=pitchFFPObj->peakLength;

	mFilterDbArr2=pitchFFPObj->mFastDbArr2;
	mFilterFreArr2=pitchFFPObj->mFastFreArr2;
	mFilterHeightArr2=pitchFFPObj->mFastHeightArr2;
	mIndexArr2=pitchFFPObj->mFastIndexArr2;
	lenArr2=pitchFFPObj->fastLenArr2;

	mFilterDbArr3=pitchFFPObj->mFastDbArr3;
	mFilterFreArr3=pitchFFPObj->mFastFreArr3;
	mFilterHeightArr3=pitchFFPObj->mFastHeightArr3;
	mIndexArr3=pitchFFPObj->mFastIndexArr3;
	lenArr3=pitchFFPObj->fastLenArr3;

	maxDBArr=pitchFFPObj->maxDBArr;

	for(int i=0;i<timeLength;i++){
		int start=0;
		int _index=0;

		len3=0;
		len2=lenArr2[i];
		maxDB=maxDBArr[i];

		// -78 filter => len2->lne3
		for(int j=0;j<len2;j++){
			if(mFilterDbArr2[i*peakLength+j]>-100){
				mFilterDbArr3[i*peakLength+len3]=mFilterDbArr2[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len3]=mFilterFreArr2[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len3]=mFilterHeightArr2[i*peakLength+j];
				mIndexArr3[i*peakLength+len3]=mIndexArr2[i*peakLength+j];

				len3++;
			}
		}

		// filter two continue >15 =>len3->len2
		{
			float _db1=0;
			float _db2=0;
			float _db3=0;
			float _db4=0;

			len2=0;
			for(int j=0;j<len3;j++){
				mFilterDbArr3[i*peakLength+len2]=mFilterDbArr3[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len2]=mFilterFreArr3[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len2]=mFilterHeightArr3[i*peakLength+j];
				mIndexArr3[i*peakLength+len2]=mIndexArr3[i*peakLength+j];
				len2++;

				if(j+3<len3){
					_db1=mFilterDbArr3[i*peakLength+j];
					_db2=mFilterDbArr3[i*peakLength+j+1];
					_db3=mFilterDbArr3[i*peakLength+j+2];
					_db4=mFilterDbArr3[i*peakLength+j+3];
					if(_db1-_db2>minDB&&_db1-_db3>minDB&&
						_db4-_db2>minDB&&_db4-_db3>minDB){ // jump

						j=j+2;
					}
				}
				// else if(j+2<len3){
				// 	_db1=mFilterDbArr3[i*peakLength+j];
				// 	_db2=mFilterDbArr3[i*peakLength+j+1];
				// 	_db3=mFilterDbArr3[i*peakLength+j+2];
				// 	if(_db1-_db2>minDB&&_db1-_db3>minDB){ // last jump
						
				// 		j=j+1;
				// 	}
				// }
			}
		}

		// left -> first <15 => len2->len3
		len3=0;
		_index=__arr_maxIndex(mFilterDbArr3+i*peakLength,len2);
		for(int j=0;j<=_index;j++){
			
			if(maxDB-mFilterDbArr3[i*peakLength+j]<minDB||
				mFilterDbArr3[i*peakLength+j]>-60||
				mFilterHeightArr3[i*peakLength+j]>18||
				__arr_has(pitchFFPObj->domIndexArr, pitchFFPObj->domLength, mIndexArr3[i*peakLength+j])){

				start=j;

				mFilterDbArr3[i*peakLength+len3]=mFilterDbArr3[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len3]=mFilterFreArr3[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len3]=mFilterHeightArr3[i*peakLength+j];
				mIndexArr3[i*peakLength+len3]=mIndexArr3[i*peakLength+j];

				len3++;
			}
		}

		// median -> relative near-filter, not cur dB-filter ???
		for(int j=start+1;j<len2-1;j++){
			if(mFilterDbArr3[i*peakLength+j-1]-mFilterDbArr3[i*peakLength+j]<minDB||
				mFilterDbArr3[i*peakLength+j+1]-mFilterDbArr3[i*peakLength+j]<minDB){

				mFilterDbArr3[i*peakLength+len3]=mFilterDbArr3[i*peakLength+j];
				mFilterFreArr3[i*peakLength+len3]=mFilterFreArr3[i*peakLength+j];
				mFilterHeightArr3[i*peakLength+len3]=mFilterHeightArr3[i*peakLength+j];
				mIndexArr3[i*peakLength+len3]=mIndexArr3[i*peakLength+j];

				len3++;
			}
		}

		// end 
		if(len2>1&&start<len2-1){ // for median fre/light
			if(mFilterDbArr3[i*peakLength+len2-2]-mFilterDbArr3[i*peakLength+len2-1]<minDB||
				len2==3||len3==2){

				mFilterDbArr3[i*peakLength+len3]=mFilterDbArr3[i*peakLength+len2-1];
				mFilterFreArr3[i*peakLength+len3]=mFilterFreArr3[i*peakLength+len2-1];
				mFilterHeightArr3[i*peakLength+len3]=mFilterHeightArr3[i*peakLength+len2-1];
				mIndexArr3[i*peakLength+len3]=mIndexArr3[i*peakLength+len2-1];

				len3++;
			}
		}
		
		lenArr3[i]=len3;
	}
}

static void __pitchFFPObj_filterFastCut(PitchFFPObj pitchFFPObj){
	int timeLength=0;
	int peakLength=0; 

	float *mFastDbArr3=NULL; // dB-filter
	float *mFastFreArr3=NULL;
	float *mFastHeightArr3=NULL;
	int *mFastIndexArr3=NULL;
	int *fastLenArr3=NULL;

	float *mFastDbArr4=NULL; // cut-4
	float *mFastFreArr4=NULL;
	float *mFastHeightArr4=NULL;
	int *mFastIndexArr4=NULL;
	int *fastLenArr4=NULL;

	timeLength=pitchFFPObj->timeLength;
	peakLength=pitchFFPObj->peakLength;

	mFastDbArr3=pitchFFPObj->mFastDbArr3;
	mFastFreArr3=pitchFFPObj->mFastFreArr3;
	mFastHeightArr3=pitchFFPObj->mFastHeightArr3;
	mFastIndexArr3=pitchFFPObj->mFastIndexArr3;
	fastLenArr3=pitchFFPObj->fastLenArr3;

	mFastDbArr4=pitchFFPObj->mFastDbArr4;
	mFastFreArr4=pitchFFPObj->mFastFreArr4;
	mFastHeightArr4=pitchFFPObj->mFastHeightArr4;
	mFastIndexArr4=pitchFFPObj->mFastIndexArr4;
	fastLenArr4=pitchFFPObj->fastLenArr4;

	for(int i=0;i<timeLength;i++){
		int len3=0;
		int len4=4;

		len3=fastLenArr3[i];

		// dB desc
		__vcorrsort1(mFastDbArr3+i*peakLength, 
					mFastFreArr3+i*peakLength,
					mFastHeightArr3+i*peakLength,
					mFastIndexArr3+i*peakLength,
					len3, 1);

		for(int j=0;j<len4;j++){
			mFastDbArr4[i*peakLength+j]=mFastDbArr3[i*peakLength+j];
			mFastFreArr4[i*peakLength+j]=mFastFreArr3[i*peakLength+j];
			mFastHeightArr4[i*peakLength+j]=mFastHeightArr3[i*peakLength+j];
			mFastIndexArr4[i*peakLength+j]=mFastIndexArr3[i*peakLength+j];
		}

		// fre asc
		__vcorrsort1(mFastFreArr4+i*peakLength, 
					mFastDbArr4+i*peakLength,
					mFastHeightArr4+i*peakLength,
					mFastIndexArr4+i*peakLength,
					len4, 0);

		fastLenArr4[i]=len4;

		// fre asc 
		__vcorrsort1(mFastFreArr3+i*peakLength, 
					mFastDbArr3+i*peakLength,
					mFastHeightArr3+i*peakLength,
					mFastIndexArr3+i*peakLength,
					len3, 0);
	}
}

static void __pitchFFPObj_stft(PitchFFPObj pitchFFPObj,float *dataArr,int dataLength,float *dbArr){
	STFTObj stftObj=NULL;

	int fftLength=0;
	int slideLength=0;
	int samplate=0;

	int peakLength=0; 

	int minIndex=0; // min/maxFre
	int maxIndex=0;

	int timeLength=0;

	// stft result ->->timeLength*(maxIndex-minIndex+1)
	float *mPowerArr=NULL;
	float *mDbArr=NULL;

	// timeLength*peakLength
	float *mPeakDbArr=NULL;
	float *mPeakFreArr=NULL;
	float *mPeakHeightArr=NULL;
	int *mIndexArr=NULL;
	int *lenArr=NULL;
	int *lowFlagArr=NULL;

	float *maxDBArr=NULL;

	// cache
	float *mRealArr=NULL; // timeLength*fftLength
	float *mImageArr=NULL;

	int len=0;
	int fLen=0;
	int rLen=0;

	float pre=0;
	float cur=0;
	float nex=0;

	float left=0;
	float right=0;

	float v1=0;
	float v2=0;

	float maxDB=0;
	float minHeight=15;

	float cutDB=-54; // -54
	float cutDB2=-58; // -60

	WindowType winType=Window_Rect;
	void (*func)(float cur,float left,float right,float *det,float *value);

	stftObj=pitchFFPObj->stftObj;

	fftLength=pitchFFPObj->fftLength;
	slideLength=pitchFFPObj->slideLength;
	samplate=pitchFFPObj->samplate;

	peakLength=pitchFFPObj->peakLength;

	minIndex=pitchFFPObj->minIndex;
	maxIndex=pitchFFPObj->maxIndex;

	timeLength=pitchFFPObj->timeLength;

	mPowerArr=pitchFFPObj->mPowerArr;
	mDbArr=pitchFFPObj->mDbArr;

	mPeakDbArr=pitchFFPObj->mPeakDbArr;
	mPeakFreArr=pitchFFPObj->mPeakFreArr;
	mPeakHeightArr=pitchFFPObj->mPeakHeightArr;
	mIndexArr=pitchFFPObj->mIndexArr;
	lenArr=pitchFFPObj->lenArr;
	lowFlagArr=pitchFFPObj->lowFlagArr;

	maxDBArr=pitchFFPObj->maxDBArr;

	mRealArr=pitchFFPObj->mRealArr;
	mImageArr=pitchFFPObj->mImageArr;
	
	winType=pitchFFPObj->winType;
	if(winType==Window_Hann){
		func=correct_hann;
	}
	else if(winType==Window_Hamm){
		func=correct_hamm;
	}
	else{ // rect
		func=correct_rect;
	}

	// 1. stft
	stftObj_stft(stftObj, dataArr, dataLength, mRealArr, mImageArr);

	// 2. power&dB
	rLen=maxIndex-minIndex+1;
	for(int i=0;i<timeLength;i++){
		for(int j=minIndex,k=0;j<=maxIndex;j++,k++){
			v1=mRealArr[i*fftLength+j];
			v2=mImageArr[i*fftLength+j];
			cur=v1*v1+v2*v2;

			mPowerArr[i*rLen+k]=cur;
			mDbArr[i*rLen+k]=10*log10f(cur/fftLength/fftLength);
		}
	}

	// 3. peak
	for(int i=0;i<timeLength;i++){
		len=0;
		fLen=0;
		maxDB=-1000;

		for(int j=1;j<rLen-1;j++){
			pre=mPowerArr[i*rLen+j-1];
			cur=mPowerArr[i*rLen+j];
			nex=mPowerArr[i*rLen+j+1];

			if(cur>pre&&cur>nex){ // peak
				float scale=0;

				float _fre=0;
				float _db=0;

				float _cur=0;
				float _pre=0;
				float _nex=0;

				float _h1=0;
				float _h2=0;
				float _height=0;

				int _index=0;

				int xFlag=0;
				int eFlag=0;
				int eFlag2=0;

				float _left=0;
				float _right=0;

				_cur=sqrtf(cur);
				_pre=sqrtf(pre);
				_nex=sqrtf(nex);

				_index=j+1;

				func(_cur,_pre,_nex,&scale,NULL);

				/***
					_db ->j is matrix
					_fre ->j+minIndex is fft'bin 
				****/
				_fre=(j+minIndex+scale)/fftLength*samplate;
				_db=mDbArr[i*rLen+j];

				// dB
				pre=mDbArr[i*rLen+j-1];
				cur=mDbArr[i*rLen+j];
				nex=mDbArr[i*rLen+j+1];

				// left height
				left=pre;
				_left=left;
				if(j-2>=0){
					left=mDbArr[i*rLen+j-2];
					_left=left;

					if(left<pre||(left>pre&&left<cur&&left-pre<2&&cur>cutDB)){
						if(j-3>=0){
							pre=mDbArr[i*rLen+j-3];
							if(pre<left){
								left=pre;
								_left=left;

								if(mDbArr[i*rLen+j-2]>mDbArr[i*rLen+j-1]&&
									mDbArr[i*rLen+j-2]<cur&&
									mDbArr[i*rLen+j-2]-mDbArr[i*rLen+j-1]<2){

									xFlag=1;
								}

								if(j-4>=0&&_db-left<minHeight&&cur>cutDB2){
									if(mDbArr[i*rLen+j-4]<pre){
										left=mDbArr[i*rLen+j-4];
										eFlag=1;
									}
								}
							}
						}
					}
					else{
						left=pre;
						_left=left;
					}
				}

				// right
				right=nex;
				_right=right;
				if(j+2<rLen){
					right=mDbArr[i*rLen+j+2];
					_right=right;

					if(right<nex||(right>nex&&right<cur&&right-nex<2&&cur>cutDB)){
						if(j+3<rLen){
							nex=mDbArr[i*rLen+j+3];

							if(nex<right){
								right=nex;
								_right=right;
								_index=j+3;

								if(j+4<rLen&&_db-right<minHeight&&!eFlag&&cur>cutDB2){
									if(mDbArr[i*rLen+j+4]<nex){
										right=mDbArr[i*rLen+j+4];
										_index=j+4;
										eFlag2=1;
									}
								}
							}
							else{
								_index=j+2;
							}
						}
					}
					else{
						right=nex;
						_right=right;
						_index=j+1;
					}
				}

				_h1=_db-left;
				_h2=_db-right;

				_height=(_h1<_h2?_h1:_h2);

				if(_height>minHeight&&xFlag&&_h1<_h2){
					mPeakDbArr[i*peakLength+len-1]=_db;
					mPeakFreArr[i*peakLength+len-1]=_fre;
					mPeakHeightArr[i*peakLength+len-1]=_height;
					mIndexArr[i*peakLength+len-1]=j;
				}
				else{
					mPeakDbArr[i*peakLength+len]=_db;
					mPeakFreArr[i*peakLength+len]=_fre;
					
					if((eFlag||eFlag2)&&cur<cutDB&&
						_height<minHeight+3){
						
						_h1=_db-_left;
						_h2=_db-_right;

						_height=(_h1<_h2?_h1:_h2);
					}
					
					mPeakHeightArr[i*peakLength+len]=_height;
					mIndexArr[i*peakLength+len]=j;

					len++;
				}

				// mPeakDbArr[i*peakLength+len]=_db;
				// mPeakFreArr[i*peakLength+len]=_fre;
				// mPeakHeightArr[i*peakLength+len]=_height;
				// mIndexArr[i*peakLength+len]=j;

				// len++;

				j=_index; // update j
			}
		}

		lenArr[i]=len;
		lowFlagArr[i]=__isLowFre(mPeakFreArr+i*peakLength, mPeakHeightArr+i*peakLength,mIndexArr+i*peakLength, len);

		// dB desc
		__vcorrsort1(mPeakDbArr+i*peakLength, 
					mPeakFreArr+i*peakLength,
					mPeakHeightArr+i*peakLength,
					mIndexArr+i*peakLength,
					len, 1);
		
		// rectify
		len=__arr_rectify(mPeakDbArr+i*peakLength,
						mPeakFreArr+i*peakLength,
						mPeakHeightArr+i*peakLength,
						mIndexArr+i*peakLength,
						len);
		lenArr[i]=len;

		maxDBArr[i]=mPeakDbArr[i*peakLength];
		if(dbArr){
			dbArr[i]=mPeakDbArr[i*peakLength];
		}
	}
}

static void __pitchFFPObj_temporal(PitchFFPObj pitchFFPObj,float *dataArr,int dataLength){
	int fftLength=0;
	int slideLength=0;

	int timeLength=0;
	float *lightArr=NULL;

	float *avgTempArr=NULL;
	float *maxTempArr=NULL;
	float *percentTempArr=NULL;

	fftLength=pitchFFPObj->fftLength;
	slideLength=pitchFFPObj->slideLength;

	timeLength=pitchFFPObj->timeLength;
	lightArr=pitchFFPObj->lightArr;

	avgTempArr=pitchFFPObj->avgTempArr;
	maxTempArr=pitchFFPObj->maxTempArr;
	percentTempArr=pitchFFPObj->percentTempArr;

	for(int i=0;i<timeLength;i++){
		lightArr[i]=__isLight(dataArr+i*slideLength, fftLength);
		maxTempArr[i]=__temproal(dataArr+i*slideLength, fftLength, pitchFFPObj->tempBase, avgTempArr+i, percentTempArr+i);
	}
}

static void __pitchFFPObj_dealData(PitchFFPObj pitchFFPObj,int dataLength){
	int fftLength=0;
	int peakLength=0;

	int minIndex=0; // min/maxFre
	int maxIndex=0;

	int timeLen=0;
	int bLen=0;

	fftLength=pitchFFPObj->fftLength;
	peakLength=pitchFFPObj->peakLength;

	minIndex=pitchFFPObj->minIndex;
	maxIndex=pitchFFPObj->maxIndex;

	bLen=maxIndex-minIndex+1;
	timeLen=stftObj_calTimeLength(pitchFFPObj->stftObj, dataLength);
	if(pitchFFPObj->timeLength<timeLen||
			pitchFFPObj->timeLength>timeLen*2){ 

		free(pitchFFPObj->mPowerArr);
		free(pitchFFPObj->mDbArr);

		free(pitchFFPObj->mPeakDbArr);
		free(pitchFFPObj->mPeakFreArr);
		free(pitchFFPObj->mPeakHeightArr);
		free(pitchFFPObj->mIndexArr);
		free(pitchFFPObj->lenArr);
		free(pitchFFPObj->lowFlagArr);

		free(pitchFFPObj->mFilterDbArr1);
		free(pitchFFPObj->mFilterFreArr1);
		free(pitchFFPObj->mFilterHeightArr1);
		free(pitchFFPObj->mIndexArr1);
		free(pitchFFPObj->lenArr1);

		free(pitchFFPObj->mFilterDbArr2);
		free(pitchFFPObj->mFilterFreArr2);
		free(pitchFFPObj->mFilterHeightArr2);
		free(pitchFFPObj->mIndexArr2);
		free(pitchFFPObj->lenArr2);

		free(pitchFFPObj->mFilterDbArr3);
		free(pitchFFPObj->mFilterFreArr3);
		free(pitchFFPObj->mFilterHeightArr3);
		free(pitchFFPObj->mIndexArr3);
		free(pitchFFPObj->lenArr3);

		free(pitchFFPObj->maxDBArr);
		free(pitchFFPObj->sucessTypeArr);
		free(pitchFFPObj->lightArr);

		free(pitchFFPObj->avgTempArr);
		free(pitchFFPObj->maxTempArr);
		free(pitchFFPObj->percentTempArr);

		free(pitchFFPObj->formatFlagArr);

		free(pitchFFPObj->formatFreArr1);
		free(pitchFFPObj->formatFreArr2);
		free(pitchFFPObj->formatFreArr3);

		free(pitchFFPObj->formatDbArr1);
		free(pitchFFPObj->formatDbArr2);
		free(pitchFFPObj->formatDbArr3);

		free(pitchFFPObj->mFastDbArr2);
		free(pitchFFPObj->mFastFreArr2);
		free(pitchFFPObj->mFastHeightArr2);
		free(pitchFFPObj->mFastIndexArr2);
		free(pitchFFPObj->fastLenArr2);

		free(pitchFFPObj->mFastDbArr3);
		free(pitchFFPObj->mFastFreArr3);
		free(pitchFFPObj->mFastHeightArr3);
		free(pitchFFPObj->mFastIndexArr3);
		free(pitchFFPObj->fastLenArr3);

		free(pitchFFPObj->mFastDbArr4);
		free(pitchFFPObj->mFastFreArr4);
		free(pitchFFPObj->mFastHeightArr4);
		free(pitchFFPObj->mFastIndexArr4);
		free(pitchFFPObj->fastLenArr4);

		pitchFFPObj->mPowerArr=__vnew(timeLen*bLen, NULL);
		pitchFFPObj->mDbArr=__vnew(timeLen*bLen, NULL);

		pitchFFPObj->mPeakDbArr=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mPeakFreArr=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mPeakHeightArr=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mIndexArr=__vnewi(timeLen*peakLength, NULL);
		pitchFFPObj->lenArr=__vnewi(timeLen, NULL);
		pitchFFPObj->lowFlagArr=__vnewi(timeLen, NULL);

		pitchFFPObj->mFilterDbArr1=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFilterFreArr1=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFilterHeightArr1=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mIndexArr1=__vnewi(timeLen*peakLength, NULL);
		pitchFFPObj->lenArr1=__vnewi(timeLen, NULL);

		pitchFFPObj->mFilterDbArr2=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFilterFreArr2=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFilterHeightArr2=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mIndexArr2=__vnewi(timeLen*peakLength, NULL);
		pitchFFPObj->lenArr2=__vnewi(timeLen, NULL);

		pitchFFPObj->mFilterDbArr3=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFilterFreArr3=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFilterHeightArr3=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mIndexArr3=__vnewi(timeLen*peakLength, NULL);
		pitchFFPObj->lenArr3=__vnewi(timeLen, NULL);

		pitchFFPObj->maxDBArr=__vnew(timeLen, NULL);
		pitchFFPObj->sucessTypeArr=__vnewi(timeLen, NULL);
		pitchFFPObj->lightArr=__vnew(timeLen, NULL);

		pitchFFPObj->avgTempArr=__vnew(timeLen, NULL);
		pitchFFPObj->maxTempArr=__vnew(timeLen, NULL);
		pitchFFPObj->percentTempArr=__vnew(timeLen, NULL);

		pitchFFPObj->formatFlagArr=__vnewi(timeLen, NULL);

		pitchFFPObj->formatFreArr1=__vnew(timeLen, NULL);
		pitchFFPObj->formatFreArr2=__vnew(timeLen, NULL);
		pitchFFPObj->formatFreArr3=__vnew(timeLen, NULL);

		pitchFFPObj->formatDbArr1=__vnew(timeLen, NULL);
		pitchFFPObj->formatDbArr2=__vnew(timeLen, NULL);
		pitchFFPObj->formatDbArr3=__vnew(timeLen, NULL);

		pitchFFPObj->mFastDbArr2=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastFreArr2=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastHeightArr2=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastIndexArr2=__vnewi(timeLen*peakLength, NULL);
		pitchFFPObj->fastLenArr2=__vnewi(timeLen, NULL);

		pitchFFPObj->mFastDbArr3=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastFreArr3=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastHeightArr3=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastIndexArr3=__vnewi(timeLen*peakLength, NULL);
		pitchFFPObj->fastLenArr3=__vnewi(timeLen, NULL);

		pitchFFPObj->mFastDbArr4=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastFreArr4=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastHeightArr4=__vnew(timeLen*peakLength, NULL);
		pitchFFPObj->mFastIndexArr4=__vnewi(timeLen*peakLength, NULL);
		pitchFFPObj->fastLenArr4=__vnewi(timeLen, NULL);

		pitchFFPObj->mRealArr=__vnew(timeLen*fftLength, NULL);
		pitchFFPObj->mImageArr=__vnew(timeLen*fftLength, NULL);
	}
	
	pitchFFPObj->timeLength=timeLen;
}

void pitchFFPObj_enableDebug(PitchFFPObj pitchFFPObj,int isDebug){

	pitchFFPObj->isDebug=isDebug;
}

void pitchFFPObj_free(PitchFFPObj pitchFFPObj){

	if(pitchFFPObj){

		stftObj_free(pitchFFPObj->stftObj);

		free(pitchFFPObj->mPowerArr);
		free(pitchFFPObj->mDbArr);

		free(pitchFFPObj->mPeakDbArr);
		free(pitchFFPObj->mPeakFreArr);
		free(pitchFFPObj->mPeakHeightArr);
		free(pitchFFPObj->mIndexArr);
		free(pitchFFPObj->lenArr);
		free(pitchFFPObj->lowFlagArr);

		free(pitchFFPObj->mFilterDbArr1);
		free(pitchFFPObj->mFilterFreArr1);
		free(pitchFFPObj->mFilterHeightArr1);
		free(pitchFFPObj->mIndexArr1);
		free(pitchFFPObj->lenArr1);

		free(pitchFFPObj->mFilterDbArr2);
		free(pitchFFPObj->mFilterFreArr2);
		free(pitchFFPObj->mFilterHeightArr2);
		free(pitchFFPObj->mIndexArr2);
		free(pitchFFPObj->lenArr2);

		free(pitchFFPObj->mFilterDbArr3);
		free(pitchFFPObj->mFilterFreArr3);
		free(pitchFFPObj->mFilterHeightArr3);
		free(pitchFFPObj->mIndexArr3);
		free(pitchFFPObj->lenArr3);

		free(pitchFFPObj->maxDBArr);
		free(pitchFFPObj->sucessTypeArr);
		free(pitchFFPObj->lightArr);

		free(pitchFFPObj->avgTempArr);
		free(pitchFFPObj->maxTempArr);
		free(pitchFFPObj->percentTempArr);

		free(pitchFFPObj->formatFlagArr);

		free(pitchFFPObj->formatFreArr1);
		free(pitchFFPObj->formatFreArr2);
		free(pitchFFPObj->formatFreArr3);

		free(pitchFFPObj->formatDbArr1);
		free(pitchFFPObj->formatDbArr2);
		free(pitchFFPObj->formatDbArr3);

		free(pitchFFPObj->mFastDbArr2);
		free(pitchFFPObj->mFastFreArr2);
		free(pitchFFPObj->mFastHeightArr2);
		free(pitchFFPObj->mFastIndexArr2);
		free(pitchFFPObj->fastLenArr2);

		free(pitchFFPObj->mFastDbArr3);
		free(pitchFFPObj->mFastFreArr3);
		free(pitchFFPObj->mFastHeightArr3);
		free(pitchFFPObj->mFastIndexArr3);
		free(pitchFFPObj->fastLenArr3);

		free(pitchFFPObj->mFastDbArr4);
		free(pitchFFPObj->mFastFreArr4);
		free(pitchFFPObj->mFastHeightArr4);
		free(pitchFFPObj->mFastIndexArr4);
		free(pitchFFPObj->fastLenArr4);

		free(pitchFFPObj->mRealArr);
		free(pitchFFPObj->mImageArr);

		free(pitchFFPObj->domIndexArr);

		free(pitchFFPObj);
	}
}

static int __isLowFre(float *freArr,float *heightArr,int *indexArr,int length){
	int flag=0;

	int num=0;

	for(int i=0;i<length-1;i++){
		int index1=0;
		int index2=0;

		if(freArr[i]<600){ // &&freArr[i+1]-freArr[i]>30 ???
			if(heightArr[i]>15&&heightArr[i+1]>15&&
				freArr[i+1]-freArr[i]>30){

				index1=indexArr[i];
				index2=indexArr[i+1];

				/***
					8seg ->62.5
					7seg ->54.6875
					6seg ->46.875
					5seg ->39.0625
				****/
				if(index2-index1<8){
					num++;
				}
			}
		}
		else{
			break;
		}
	}

	if(num>=4){
		flag=1;
	}

	return flag;
}

/***
	max<-18 & >0.4
****/
static float __isLight(float *dataArr,int length){
	float result=0;

	int count=0;
	float value=0;
	
	if(!dataArr||!length){
		return 0;
	}

	for(int i=0;i<length;i++){
		value=20*log10f(fabsf(dataArr[i])+1e-8);
		if(value>-18){
			return 0;
		}

		if(value>-24){
			count++;
		}
	}

	result=1.0*(length-count)/length;

	return result;
}

static float __temproal(float *dataArr,int length,float base,float *avgValue,float *percentValue){
	float max=-100;

	int count=0;
	float value=0;
	float total=0;
	
	if(!dataArr||!length){
		return 0;
	}

	for(int i=0;i<length;i++){
		value=20*log10f(fabsf(dataArr[i])+1e-8);
		if(value<-36){
			value=-36;
		}

		if(max<value){
			max=value;
		}

		total+=value;

		if(value>-base){
			count++;
		}
	}

	*avgValue=total/length;
	*percentValue=1.0*(length-count)/length;

	return max;
}

static int __arr_rectify(float *dbArr,float *freArr,float *heightArr,int *indexArr,int length){
	int len=0;
	int offset=0;

	if(length<3){
		return length;
	}

	len=length;
	if(abs(indexArr[0]-indexArr[1])<=4&&
		dbArr[0]-dbArr[1]<3){

		float s1=0;
		float s2=0;

		s1=fabsf(2*freArr[0]-freArr[2]);
		s2=fabsf(2*freArr[1]-freArr[2]);
		if(s1<s2){
			offset=1;
		}
		else{
			offset=0;
		}

		for(int i=offset;i<length-1;i++){
			dbArr[i]=dbArr[i+1];
			freArr[i]=freArr[i+1];
			heightArr[i]=heightArr[i+1];
			indexArr[i]=indexArr[i+1];
		}

		len=length-1;
	}
	else if(abs(indexArr[1]-indexArr[2])<=4&&
		dbArr[1]-dbArr[2]<3){

		float s1=0;
		float s2=0;

		if(freArr[0]>freArr[1]){
			s1=fabsf(2*freArr[1]-freArr[0]);
			s2=fabsf(2*freArr[2]-freArr[0]);
		}
		else{
			s1=fabsf(freArr[1]-2*freArr[0]);
			s2=fabsf(freArr[2]-2*freArr[0]);
		}

		if(s1<s2){
			offset=2;
		}
		else{
			offset=1;
		}

		for(int i=offset;i<length-1;i++){
			dbArr[i]=dbArr[i+1];
			freArr[i]=freArr[i+1];
			heightArr[i]=heightArr[i+1];
			indexArr[i]=indexArr[i+1];
		}

		len=length-1;
	}
	else if(abs(indexArr[0]-indexArr[2])<=4&&
		dbArr[0]-dbArr[2]<3){

		float s1=0;
		float s2=0;

		s1=fabsf(2*freArr[0]-freArr[1]);
		s2=fabsf(2*freArr[2]-freArr[1]);
		if(s1<s2){
			offset=2;
		}
		else{
			offset=0;
		}

		for(int i=offset;i<length-1;i++){
			dbArr[i]=dbArr[i+1];
			freArr[i]=freArr[i+1];
			heightArr[i]=heightArr[i+1];
			indexArr[i]=indexArr[i+1];
		}

		len=length-1;
	}

	return len;
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

static int __arr_has(int *arr,int length,int value){
	int flag=0;

	for(int i=0;i<length;i++){
		if(arr[i]==value){
			flag=1;
			break;
		}
	}

	return flag;
}










