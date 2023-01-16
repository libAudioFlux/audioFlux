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

#include "../classic/trist.h"

#include "../stft_algorithm.h"

#include "_pitch_stft.h"

struct OpaquePitchSTFT{
	int isContinue;

	STFTObj stftObj;

	int fftLength;
	int slideLength;

	int cutLength; // 12
	int peakLength; // (maxIndex-minIndex+1)/2+1

	int minIndex; // min/maxFre
	int maxIndex;

	int timeLength;

	// stft result ->timeLength*peakLength
	float *mDbArr; 
	float *mCorrectFreArr; 
	float *mHeightArr;
	int *mMidiArr; // tone&correctFre

	float *mFeatureFreArr; // count1+count2
	float *mFeatureDbArr;
	float *mFeatureHeightArr;
	int *mFeatureMidiArr;

	int *countArr1; // 1000 height>=20 count
	int *countArr2;	

	int *lenArr;

	// cache
	float *mRealArr; // timeLength*fftLength
	float *mImageArr;

	int samplate;
	WindowType winType;

	int isDebug;
};

static void __pitchSTFTObj_dealData(PitchSTFTObj pitchSTFTObj,int dataLength);

static void __pitchSTFTObj_stft(PitchSTFTObj pitchSTFTObj,float *dataArr,int dataLength,float *dbArr);
static void __pitchSTFTObj_sub(PitchSTFTObj pitchSTFTObj,float *freArr);

static int __trist_sub(float *corArr,float *dbArr,float *heightArr,int *midiArr,int length,float *outFre);
// 0 desc 1 asc
static void __arr_relateSort(float *arr1,float *arr2,float *arr3,int *arr4,int length,int type);

static float __hanScaleCorrect(float cur,float left,float right);

/***
	samplate 32000
	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int pitchSTFTObj_new(PitchSTFTObj *pitchSTFTObj,
				int *samplate,float *lowFre,float *highFre,
				int *radix2Exp,int *slideLength,WindowType *windowType,
				int *isContinue){
	int status=0;

	int _samplate=32000;
	float _lowFre=27; // 27.5
	float _highFre=2094; // 2093
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
	PitchSTFTObj pitch=NULL;

	pitch=*pitchSTFTObj=(PitchSTFTObj )calloc(1,sizeof(struct OpaquePitchSTFT ));

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

	minIndex=floorf(_highFre*fftLength/_samplate);
	maxIndex=ceilf(_lowFre*fftLength/_samplate);
	if(maxIndex>=fftLength/2){
		maxIndex=fftLength/2-1;
	}

	if(minIndex>=maxIndex){
		minIndex=3;
		maxIndex=ceilf(2093*fftLength/_samplate);
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
	
	return status;
}

int pitchSTFTObj_calTimeLength(PitchSTFTObj pitchSTFTObj,int dataLength){
	int timeLen=0;

	timeLen=stftObj_calTimeLength(pitchSTFTObj->stftObj, dataLength);
	return timeLen;
}

void pitchSTFTObj_pitch(PitchSTFTObj pitchSTFTObj,float *dataArr,int dataLength,
					float *freArr,float *dbArr){

	// 1. dealData
	__pitchSTFTObj_dealData(pitchSTFTObj,dataLength);
	// 2. stft
	__pitchSTFTObj_stft(pitchSTFTObj, dataArr, dataLength,dbArr);
	// 3. sub
	__pitchSTFTObj_sub(pitchSTFTObj,freArr);

}

int pitchSTFTObj_getCorrData(PitchSTFTObj pitchSTFTObj,float **mCorrArr,int **lenArr){
	int mLen=0;

	mLen=pitchSTFTObj->peakLength;
	if(mCorrArr){
		*mCorrArr=pitchSTFTObj->mCorrectFreArr;
	}

	if(lenArr){
		*lenArr=pitchSTFTObj->lenArr;
	}

	return mLen;
}

static void __pitchSTFTObj_sub(PitchSTFTObj pitchSTFTObj,float *freArr){
	float *mDbArr=NULL; 
	float *mCorrectFreArr=NULL; 
	float *mHeightArr=NULL; 
	int *mMidiArr=NULL;	

	float *mFeatureFreArr=NULL;
	float *mFeatureDbArr=NULL;
	float *mFeatureHeightArr=NULL;
	int *mFeatureMidiArr=NULL;

	int *countArr1=NULL;
	int *countArr2=NULL;

	int *lenArr=NULL;

	int cutLength=0; // 12
	int peakLength=0;

	int timeLength=0;
	int len=0;

	float *corrArr=NULL;
	float *dbArr1=NULL;
	float *heightArr1=NULL;
	int *midiArr1=NULL;

	float *featureArr=NULL;
	float *dbArr2=NULL;
	float *heightArr2=NULL;
	int *midiArr2=NULL;

	int flag=0;

	mDbArr=pitchSTFTObj->mDbArr;
	mCorrectFreArr=pitchSTFTObj->mCorrectFreArr;
	mHeightArr=pitchSTFTObj->mHeightArr;
	mMidiArr=pitchSTFTObj->mMidiArr;

	mFeatureFreArr=pitchSTFTObj->mFeatureFreArr;
	mFeatureDbArr=pitchSTFTObj->mFeatureDbArr;
	mFeatureHeightArr=pitchSTFTObj->mFeatureHeightArr;
	mFeatureMidiArr=pitchSTFTObj->mFeatureMidiArr;

	countArr1=pitchSTFTObj->countArr1;
	countArr2=pitchSTFTObj->countArr2;

	lenArr=pitchSTFTObj->lenArr;

	cutLength=pitchSTFTObj->cutLength;
	peakLength=pitchSTFTObj->peakLength;

	timeLength=pitchSTFTObj->timeLength;
	for(int i=0;i<timeLength;i++){
		len=lenArr[i];

		corrArr=mCorrectFreArr+i*peakLength;
		dbArr1=mDbArr+i*peakLength;
		heightArr1=mHeightArr+i*peakLength;
		midiArr1=mMidiArr+i*peakLength;

		featureArr=mFeatureFreArr+i*peakLength;
		dbArr2=mFeatureDbArr+i*peakLength;
		heightArr2=mFeatureHeightArr+i*peakLength;
		midiArr2=mFeatureMidiArr+i*peakLength;
		
		flag=trist(corrArr,dbArr1,heightArr1,midiArr1,len,
				featureArr,dbArr2,heightArr2,midiArr2,countArr1[i],countArr2[i],
				freArr+i);
	}
}

static void __pitchSTFTObj_stft(PitchSTFTObj pitchSTFTObj,float *dataArr,int dataLength,float *dbArr){
	STFTObj stftObj=NULL;

	int fftLength=0;
	int slideLength=0;
	int samplate=0;

	int peakLength=0; 

	int minIndex=0; // min/maxFre
	int maxIndex=0;

	int timeLength=0;

	// stft result ->timeLength*cutLength
	float *mDbArr=NULL; 
	float *mCorrectFreArr=NULL; 
	float *mHeightArr=NULL;
	int *mMidiArr=NULL;

	float *mFeatureFreArr=NULL; // count1+count2
	float *mFeatureDbArr=NULL;
	float *mFeatureHeightArr=NULL;
	int *mFeatureMidiArr=NULL;

	int *countArr1=NULL; // 1000 height>=20 count
	int *countArr2=NULL;	
	int *lenArr=NULL;

	float minHeight=20;
	float compareValue=2.6;

	// cache
	float *mRealArr=NULL; // timeLength*fftLength
	float *mImageArr=NULL;

	int len=0;
	int fLen=0;

	float pre=0;
	float cur=0;
	float nex=0;

	float left=0;
	float right=0;

	float v1=0;
	float v2=0;

	int sIndex=0;

	WindowType winType=Window_Rect;
	void (*func)(float cur,float left,float right,float *det,float *value);

	stftObj=pitchSTFTObj->stftObj;

	fftLength=pitchSTFTObj->fftLength;
	slideLength=pitchSTFTObj->slideLength;
	samplate=pitchSTFTObj->samplate;

	peakLength=pitchSTFTObj->peakLength;

	minIndex=pitchSTFTObj->minIndex;
	maxIndex=pitchSTFTObj->maxIndex;

	timeLength=pitchSTFTObj->timeLength;

	mDbArr=pitchSTFTObj->mDbArr;
	mCorrectFreArr=pitchSTFTObj->mCorrectFreArr;
	mHeightArr=pitchSTFTObj->mHeightArr;
	mMidiArr=pitchSTFTObj->mMidiArr;

	mFeatureFreArr=pitchSTFTObj->mFeatureFreArr;
	mFeatureDbArr=pitchSTFTObj->mFeatureDbArr;
	mFeatureHeightArr=pitchSTFTObj->mFeatureHeightArr;
	mFeatureMidiArr=pitchSTFTObj->mFeatureMidiArr;

	countArr1=pitchSTFTObj->countArr1;
	countArr2=pitchSTFTObj->countArr2;
	lenArr=pitchSTFTObj->lenArr;

	mRealArr=pitchSTFTObj->mRealArr;
	mImageArr=pitchSTFTObj->mImageArr;
	
	winType=pitchSTFTObj->winType;
	if(winType==Window_Hann){
		func=correct_hann;
	}
	else if(winType==Window_Hamm){
		func=correct_hamm;
	}
	else{ // rect
		func=correct_rect;
	}

	sIndex=roundf(1000.0*fftLength/samplate);

	// 1. stft
	stftObj_stft(stftObj, dataArr, dataLength, mRealArr, mImageArr);

	// 2. dB/fre
	for(int i=0;i<timeLength;i++){
		len=0;
		fLen=0;
		countArr1[i]=0;
		countArr2[i]=0;
		for(int j=minIndex+1;j<maxIndex;j++){
			v1=mRealArr[i*fftLength+j-1];
			v2=mImageArr[i*fftLength+j-1];
			pre=v1*v1+v2*v2;

			v1=mRealArr[i*fftLength+j];
			v2=mImageArr[i*fftLength+j];
			cur=v1*v1+v2*v2;

			v1=mRealArr[i*fftLength+j+1];
			v2=mImageArr[i*fftLength+j+1];
			nex=v1*v1+v2*v2;

			if(cur>pre&&cur>nex){ // peak
				float scale=0;

				float _fre=0;
				float _db=0;
				int _midi=0;

				float _cur=0;
				float _pre=0;
				float _nex=0;

				float _h1=0;
				float _h2=0;
				float _height=0;

				int _index=0;

				_cur=sqrtf(cur);
				_pre=sqrtf(pre);
				_nex=sqrtf(nex);

				_index=j+1;

				func(_cur,_pre,_nex,&scale,NULL);
				// printf("scale1 is %.3f\n",scale);
				// scale=__hanScaleCorrect(_cur,_pre,_nex);
				// printf("scale2 is %.3f\n",scale);

				_fre=(j+scale)/fftLength*samplate;
				_db=10*log10f(cur/fftLength/fftLength);
				_midi=util_freToMidi(_fre);

				mCorrectFreArr[i*peakLength+len]=_fre;
				mDbArr[i*peakLength+len]=_db;
				mMidiArr[i*peakLength+len]=_midi;

				// left height
				left=pre;
				if(j-2>=0){
					v1=mRealArr[i*fftLength+j-2];
					v2=mImageArr[i*fftLength+j-2];
					left=v1*v1+v2*v2;

					if(left<pre){
						if(j-3>=0){
							v1=mRealArr[i*fftLength+j-3];
							v2=mImageArr[i*fftLength+j-3];
							pre=v1*v1+v2*v2;

							if(pre<left){
								left=pre;
								// if(j-4>=0){
								// 	v1=mRealArr[i*fftLength+j-4];
								// 	v2=mImageArr[i*fftLength+j-4];
								// 	pre=v1*v1+v2*v2;

								// 	if(pre<left){
								// 		left=pre;
								// 	}
								// }
							}
						}
					}
					else{
						left=pre;
					}
				}

				// right
				right=nex;
				if(j+2<fftLength/2){
					v1=mRealArr[i*fftLength+j+2];
					v2=mImageArr[i*fftLength+j+2];
					right=v1*v1+v2*v2;

					if(right<nex){
						if(j+3<fftLength/2){
							v1=mRealArr[i*fftLength+j+3];
							v2=mImageArr[i*fftLength+j+3];
							nex=v1*v1+v2*v2;

							if(nex<right){
								right=nex;
								// if(j+4<fftLength/2){
								// 	v1=mRealArr[i*fftLength+j+4];
								// 	v2=mImageArr[i*fftLength+j+4];
								// 	nex=v1*v1+v2*v2;

								// 	if(nex<right){
								// 		right=nex;
								// 		_index=j+4;
								// 	}
								// }
								// else{
								// 	_index=j+3;
								// }

								_index=j+3;
							}
							else{
								_index=j+2;
							}
						}
					}
					else{
						right=nex;
						_index=j+1;
					}
				}

				_h1=_db-10*log10f(left/fftLength/fftLength);
				_h2=_db-10*log10f(right/fftLength/fftLength);

				_height=(_h1<_h2?_h1:_h2);
				mHeightArr[i*peakLength+len]=_height;

				// countArr1/countArr2
				if(_height>=minHeight){
					if(j<sIndex){
						countArr1[i]++;
					}
					else if(j<2*sIndex){
						countArr2[i]++;
					}

					mFeatureFreArr[i*peakLength+fLen]=_fre;
					mFeatureDbArr[i*peakLength+fLen]=_db;
					mFeatureHeightArr[i*peakLength+fLen]=_height;
					mFeatureMidiArr[i*peakLength+fLen]=_midi;

					fLen++;
				}

				len++;
				j=_index; // update j
			}
		}

		__arr_relateSort(mDbArr+i*peakLength, 
						mCorrectFreArr+i*peakLength,
						mHeightArr+i*peakLength,
						mMidiArr+i*peakLength,
						len, 0);
		if(dbArr){
			dbArr[i]=mDbArr[i*peakLength];
		}
		lenArr[i]=len;
	}

	// 3. debug
	if(pitchSTFTObj->isDebug){
		printf("stft start !!!\n");
		for(int i=0;i<timeLength;i++){
			printf("index[%d], %0.3f:\n",i,1.0*i*slideLength/samplate);
			int len=0;
			int fLen=0;

			len=lenArr[i];
			fLen=countArr1[i]+countArr2[i];

			printf("	freArr:\n");
			printf("		");
			for(int j=0;j<len;j++){
				printf("[%d]:%0.1f,",j,mCorrectFreArr[i*peakLength+j]);
			}
			printf("\n");

			printf("	dbArr:\n");
			printf("		");
			for(int j=0;j<len;j++){
				printf("[%d]:%0.1f,",j,mDbArr[i*peakLength+j]);
			}
			printf("\n");

			printf("	heightArr:\n");
			printf("		");
			for(int j=0;j<len;j++){
				printf("[%d]:%0.1f,",j,mHeightArr[i*peakLength+j]);
			}
			printf("\n");

			printf("	count: %d, %d",countArr1[i],countArr2[i]);
			printf("\n");
			printf("\n");

			printf("	freArr1:\n");
			printf("		");
			for(int j=0;j<fLen;j++){
				printf("[%d]:%0.1f,",j,mFeatureFreArr[i*peakLength+j]);
			}
			printf("\n");

			printf("	dbArr1:\n");
			printf("		");
			for(int j=0;j<fLen;j++){
				printf("[%d]:%0.1f,",j,mFeatureDbArr[i*peakLength+j]);
			}
			printf("\n");

			printf("	heightArr1:\n");
			printf("		");
			for(int j=0;j<fLen;j++){
				printf("[%d]:%0.1f,",j,mFeatureHeightArr[i*peakLength+j]);
			}
			printf("\n");
			printf("\n");
		}
		printf("stft end !!!\n");
	}
}

static void __pitchSTFTObj_dealData(PitchSTFTObj pitchSTFTObj,int dataLength){
	int fftLength=0;
	int peakLength=0;

	float *mDbArr=NULL; 
	float *mCorrectFreArr=NULL; 
	int *lenArr=NULL;

	float *mRealArr=NULL;
	float *mImageArr=NULL;

	int minIndex=0; // min/maxFre
	int maxIndex=0;

	int timeLen=0;

	fftLength=pitchSTFTObj->fftLength;
	peakLength=pitchSTFTObj->peakLength;

	mDbArr=pitchSTFTObj->mDbArr;
	mCorrectFreArr=pitchSTFTObj->mCorrectFreArr;
	lenArr=pitchSTFTObj->lenArr;

	mRealArr=pitchSTFTObj->mRealArr;
	mImageArr=pitchSTFTObj->mImageArr;

	minIndex=pitchSTFTObj->minIndex;
	maxIndex=pitchSTFTObj->maxIndex;

	timeLen=stftObj_calTimeLength(pitchSTFTObj->stftObj, dataLength);
	if(pitchSTFTObj->timeLength<timeLen||
			pitchSTFTObj->timeLength>timeLen*2){ 
		free(pitchSTFTObj->mDbArr);
		free(pitchSTFTObj->mCorrectFreArr);
		free(pitchSTFTObj->mHeightArr);
		free(pitchSTFTObj->mMidiArr);

		free(pitchSTFTObj->mFeatureFreArr);
		free(pitchSTFTObj->mFeatureDbArr);
		free(pitchSTFTObj->mFeatureHeightArr);
		free(pitchSTFTObj->mFeatureMidiArr);

		free(pitchSTFTObj->countArr1);
		free(pitchSTFTObj->countArr2);
		free(pitchSTFTObj->lenArr);

		pitchSTFTObj->mDbArr=__vnew(timeLen*peakLength, NULL);
		pitchSTFTObj->mCorrectFreArr=__vnew(timeLen*peakLength, NULL);
		pitchSTFTObj->mHeightArr=__vnew(timeLen*peakLength, NULL);
		pitchSTFTObj->mMidiArr=__vnewi(timeLen*peakLength, NULL);

		pitchSTFTObj->mFeatureFreArr=__vnew(timeLen*peakLength, NULL);
		pitchSTFTObj->mFeatureDbArr=__vnew(timeLen*peakLength, NULL);
		pitchSTFTObj->mFeatureHeightArr=__vnew(timeLen*peakLength, NULL);
		pitchSTFTObj->mFeatureMidiArr=__vnewi(timeLen*peakLength, NULL);

		pitchSTFTObj->countArr1=__vnewi(timeLen, NULL);
		pitchSTFTObj->countArr2=__vnewi(timeLen, NULL);
		pitchSTFTObj->lenArr=__vnewi(timeLen, NULL);

		pitchSTFTObj->mRealArr=__vnew(timeLen*fftLength, NULL);
		pitchSTFTObj->mImageArr=__vnew(timeLen*fftLength, NULL);
	}
	
	pitchSTFTObj->timeLength=timeLen;
}

void pitchSTFTObj_enableDebug(PitchSTFTObj pitchSTFTObj,int isDebug){

	pitchSTFTObj->isDebug=isDebug;
}

void pitchSTFTObj_free(PitchSTFTObj pitchSTFTObj){

	if(pitchSTFTObj){

		stftObj_free(pitchSTFTObj->stftObj);

		free(pitchSTFTObj->mDbArr);
		free(pitchSTFTObj->mCorrectFreArr);
		free(pitchSTFTObj->mHeightArr);
		free(pitchSTFTObj->mMidiArr);

		free(pitchSTFTObj->mFeatureFreArr);
		free(pitchSTFTObj->mFeatureDbArr);
		free(pitchSTFTObj->mFeatureHeightArr);
		free(pitchSTFTObj->mFeatureMidiArr);

		free(pitchSTFTObj->countArr1);
		free(pitchSTFTObj->countArr2);
		free(pitchSTFTObj->lenArr);

		free(pitchSTFTObj->mRealArr);
		free(pitchSTFTObj->mImageArr);

		free(pitchSTFTObj);
	}
}

// 0 desc 1 asc
static void __arr_relateSort(float *arr1,float *arr2,float *arr3,int *arr4,int length,int type){
	float value1=0;
	float value2=0;
	float value3=0;
	int value4=0;

	for(int i=0;i<length;i++){
		for(int j=i+1;j<length;j++){
			if(type==0){ // 降序
				if(arr1[i]<arr1[j]){
					value1=arr1[i];
					arr1[i]=arr1[j];
					arr1[j]=value1;

					if(arr2){
						value2=arr2[i];
						arr2[i]=arr2[j];
						arr2[j]=value2;
					}
					
					if(arr3){
						value3=arr3[i];
						arr3[i]=arr3[j];
						arr3[j]=value3;
					}

					if(arr4){
						value4=arr4[i];
						arr4[i]=arr4[j];
						arr4[j]=value4;
					}
				}
			}
			else{ // 升序
				if(arr1[i]>arr1[j]){
					value1=arr1[i];
					arr1[i]=arr1[j];
					arr1[j]=value1;

					if(arr2){
						value2=arr2[i];
						arr2[i]=arr2[j];
						arr2[j]=value2;
					}
					
					if(arr3){
						value3=arr3[i];
						arr3[i]=arr3[j];
						arr3[j]=value3;
					}

					if(arr4){
						value4=arr4[i];
						arr4[i]=arr4[j];
						arr4[j]=value4;
					}
				}
			}
		}
	}
}

float __hanScaleCorrect(float cur,float left,float right){
    float d1=0;
    float d2=0;
    
    d1=right-left;
    d2=(d1>=0)?(2*right-cur)/(cur+right):(cur-2*left)/(cur+left);
    
    return d2;
}








