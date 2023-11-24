// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/flux_window.h"
#include "dsp/flux_correct.h"
#include "dsp/fft_algorithm.h"

#include "stft_algorithm.h"

#include "harmonic_algorithm.h"

struct OpaqueHarmonic{
	int isContinue;

	STFTObj stftObj;

	int fftLength;
	int slideLength;

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

	// cache
	float *mRealArr; // timeLength*fftLength
	float *mImageArr;

	int samplate;
	WindowType winType;

	int isDebug;
};

static void __harmonicObj_dealData(HarmonicObj harmonicObj,int dataLength);

static void __harmonicObj_stft(HarmonicObj harmonicObj,float *dataArr,int dataLength);

static void __harmonicObj_filter(HarmonicObj harmonicObj);

static void __harmonicObj_filterHeight(HarmonicObj harmonicObj);
static void __harmonicObj_filterNear(HarmonicObj harmonicObj);
static void __harmonicObj_filterDB(HarmonicObj harmonicObj);

static int __arr_maxIndex(float *arr,int length);

int harmonicObj_new(HarmonicObj *harmonicObj,
					int *samplate,
					int *radix2Exp,WindowType *windowType,int *slideLength){
	int status=0;

	int _samplate=32000;

	int _radix2Exp=12;
	int _slideLength=0;
	WindowType _winType=Window_Hamm;

	int fftLength=0;
	int peakLength=0;

	int minIndex=0; // min/maxFre
	int maxIndex=0;

	STFTObj stftObj=NULL;
	HarmonicObj hr=NULL;

	hr=*harmonicObj=(HarmonicObj )calloc(1,sizeof(struct OpaqueHarmonic ));

	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
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

	maxIndex=fftLength/2-1;
	peakLength=(maxIndex-minIndex)/2+1;

	stftObj_new(&stftObj, _radix2Exp, &_winType, &_slideLength, NULL);

	hr->stftObj=stftObj;

	hr->fftLength=fftLength;
	hr->slideLength=_slideLength;

	hr->peakLength=peakLength;

	hr->minIndex=minIndex;
	hr->maxIndex=maxIndex;

	hr->samplate=_samplate;
	hr->winType=_winType;
	
	return status;
}

int harmonicObj_calTimeLength(HarmonicObj harmonicObj,int dataLength){
	int timeLen=0;

	timeLen=stftObj_calTimeLength(harmonicObj->stftObj, dataLength);
	return timeLen;
}

void harmonicObj_exec(HarmonicObj harmonicObj,float *dataArr,int dataLength){
	// 1. dealData
	__harmonicObj_dealData(harmonicObj,dataLength);
	// 2. stft
	__harmonicObj_stft(harmonicObj, dataArr, dataLength);
	// 3. filter
	__harmonicObj_filter(harmonicObj);
}

void harmonicObj_harmonicCount(HarmonicObj harmonicObj,float low,float high,int *countArr){
	int timeLength=0;
	int peakLength=0; 

	float *mFilterDbArr3=NULL; // dB-filter
	float *mFilterFreArr3=NULL; 
	float *mFilterHeightArr3=NULL;
	int *mIndexArr3=NULL;
	int *lenArr3=NULL;

	float *freArr=NULL;
	int len3=0;

	mFilterDbArr3=harmonicObj->mFilterDbArr3;
	mFilterFreArr3=harmonicObj->mFilterFreArr3;
	mFilterHeightArr3=harmonicObj->mFilterHeightArr3;
	mIndexArr3=harmonicObj->mIndexArr3;
	lenArr3=harmonicObj->lenArr3;

	timeLength=harmonicObj->timeLength;
	peakLength=harmonicObj->peakLength;

	for(int i=0;i<timeLength;i++){
		freArr=mFilterFreArr3+i*peakLength;
		len3=lenArr3[i];

		for(int j=0;j<len3;j++){
			if(freArr[j]>=high){
				break;
			}

			if(freArr[j]>low&&freArr[j]<high){
				countArr[i]++;
			}
		}
	}
}

static void __harmonicObj_dealData(HarmonicObj harmonicObj,int dataLength){
	int fftLength=0;
	int peakLength=0;

	int minIndex=0; // min/maxFre
	int maxIndex=0;

	int timeLen=0;
	int bLen=0;

	fftLength=harmonicObj->fftLength;
	peakLength=harmonicObj->peakLength;

	minIndex=harmonicObj->minIndex;
	maxIndex=harmonicObj->maxIndex;

	bLen=maxIndex-minIndex+1;
	timeLen=stftObj_calTimeLength(harmonicObj->stftObj, dataLength);
	if(harmonicObj->timeLength<timeLen||
			harmonicObj->timeLength>timeLen*2){ 

		free(harmonicObj->mPowerArr);
		free(harmonicObj->mDbArr);

		free(harmonicObj->mPeakDbArr);
		free(harmonicObj->mPeakFreArr);
		free(harmonicObj->mPeakHeightArr);
		free(harmonicObj->mIndexArr);
		free(harmonicObj->lenArr);

		free(harmonicObj->mFilterDbArr1);
		free(harmonicObj->mFilterFreArr1);
		free(harmonicObj->mFilterHeightArr1);
		free(harmonicObj->mIndexArr1);
		free(harmonicObj->lenArr1);

		free(harmonicObj->mFilterDbArr2);
		free(harmonicObj->mFilterFreArr2);
		free(harmonicObj->mFilterHeightArr2);
		free(harmonicObj->mIndexArr2);
		free(harmonicObj->lenArr2);

		free(harmonicObj->mFilterDbArr3);
		free(harmonicObj->mFilterFreArr3);
		free(harmonicObj->mFilterHeightArr3);
		free(harmonicObj->mIndexArr3);
		free(harmonicObj->lenArr3);

		free(harmonicObj->maxDBArr);

		harmonicObj->mPowerArr=__vnew(timeLen*bLen, NULL);
		harmonicObj->mDbArr=__vnew(timeLen*bLen, NULL);

		harmonicObj->mPeakDbArr=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mPeakFreArr=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mPeakHeightArr=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mIndexArr=__vnewi(timeLen*peakLength, NULL);
		harmonicObj->lenArr=__vnewi(timeLen, NULL);

		harmonicObj->mFilterDbArr1=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mFilterFreArr1=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mFilterHeightArr1=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mIndexArr1=__vnewi(timeLen*peakLength, NULL);
		harmonicObj->lenArr1=__vnewi(timeLen, NULL);

		harmonicObj->mFilterDbArr2=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mFilterFreArr2=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mFilterHeightArr2=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mIndexArr2=__vnewi(timeLen*peakLength, NULL);
		harmonicObj->lenArr2=__vnewi(timeLen, NULL);

		harmonicObj->mFilterDbArr3=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mFilterFreArr3=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mFilterHeightArr3=__vnew(timeLen*peakLength, NULL);
		harmonicObj->mIndexArr3=__vnewi(timeLen*peakLength, NULL);
		harmonicObj->lenArr3=__vnewi(timeLen, NULL);

		harmonicObj->maxDBArr=__vnew(timeLen, NULL);

		harmonicObj->mRealArr=__vnew(timeLen*fftLength, NULL);
		harmonicObj->mImageArr=__vnew(timeLen*fftLength, NULL);
	}
	
	harmonicObj->timeLength=timeLen;
}

static void __harmonicObj_stft(HarmonicObj harmonicObj,float *dataArr,int dataLength){
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

	float minHeight=15;
	float cutDB=-50;

	WindowType winType=Window_Rect;
	void (*func)(float cur,float left,float right,float *det,float *value);

	stftObj=harmonicObj->stftObj;

	fftLength=harmonicObj->fftLength;
	slideLength=harmonicObj->slideLength;
	samplate=harmonicObj->samplate;

	peakLength=harmonicObj->peakLength;

	minIndex=harmonicObj->minIndex;
	maxIndex=harmonicObj->maxIndex;

	timeLength=harmonicObj->timeLength;

	mPowerArr=harmonicObj->mPowerArr;
	mDbArr=harmonicObj->mDbArr;

	mPeakDbArr=harmonicObj->mPeakDbArr;
	mPeakFreArr=harmonicObj->mPeakFreArr;
	mPeakHeightArr=harmonicObj->mPeakHeightArr;
	mIndexArr=harmonicObj->mIndexArr;
	lenArr=harmonicObj->lenArr;

	maxDBArr=harmonicObj->maxDBArr;

	mRealArr=harmonicObj->mRealArr;
	mImageArr=harmonicObj->mImageArr;
	
	winType=harmonicObj->winType;
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

		for(int j=1;j<rLen-1;j++){
			pre=mPowerArr[i*rLen+j-1];
			cur=mPowerArr[i*rLen+j];
			nex=mPowerArr[i*rLen+j+1];

			if(cur>pre&&cur>nex){ // peak
				float scale=0;

				float _fre=0;
				float _db=0;

				float _h1=0;
				float _h2=0;
				float _height=0;

				int _index=0;

				int xFlag=0;
				int eFlag=0;

				// float _cur=0;
				// float _pre=0;
				// float _nex=0;

				// _cur=sqrtf(cur);
				// _pre=sqrtf(pre);
				// _nex=sqrtf(nex);

				// func(_cur,_pre,_nex,&scale,NULL);

				_index=j+1;

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
				if(j-2>=0){
					left=mDbArr[i*rLen+j-2];

					if(left<pre||(left>pre&&left<cur&&left-pre<2&&cur>cutDB)){
						if(j-3>=0){
							pre=mDbArr[i*rLen+j-3];
							if(pre<left){
								left=pre;

								if(mDbArr[i*rLen+j-2]>mDbArr[i*rLen+j-1]&&
									mDbArr[i*rLen+j-2]<cur&&
									mDbArr[i*rLen+j-2]-mDbArr[i*rLen+j-1]<2){

									xFlag=1;
								}

								if(j-4>=0&&_db-left<minHeight&&cur>cutDB){
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
					}
				}

				// right
				right=nex;
				if(j+2<rLen){
					right=mDbArr[i*rLen+j+2];;

					if(right<nex||(right>nex&&right<cur&&right-nex<2&&cur>cutDB)){
						if(j+3<rLen){
							nex=mDbArr[i*rLen+j+3];

							if(nex<right){
								right=nex;
								_index=j+3;

								if(j+4<rLen&&_db-right<minHeight&&!eFlag&&cur>cutDB){
									if(mDbArr[i*rLen+j+4]<nex){
										right=mDbArr[i*rLen+j+4];
										_index=j+4;
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
					mPeakHeightArr[i*peakLength+len]=_height;
					mIndexArr[i*peakLength+len]=j;

					len++;
				}

				j=_index; // update j
			}
		}

		lenArr[i]=len;

		// dB desc
		__vcorrsort1(mPeakDbArr+i*peakLength, 
					mPeakFreArr+i*peakLength,
					mPeakHeightArr+i*peakLength,
					mIndexArr+i*peakLength,
					len, 1);

		maxDBArr[i]=mPeakDbArr[i*peakLength];
	}
}

static void __harmonicObj_filter(HarmonicObj harmonicObj){

	__harmonicObj_filterHeight(harmonicObj);
	__harmonicObj_filterNear(harmonicObj);
	__harmonicObj_filterDB(harmonicObj);

}

static void __harmonicObj_filterHeight(HarmonicObj harmonicObj){
	int timeLength=0;
	int peakLength=0; 

	float *mPeakDbArr=NULL;
	float *mPeakFreArr=NULL;
	float *mPeakHeightArr=NULL;
	int *mIndexArr=NULL;
	int *lenArr=NULL;

	float *mFilterDbArr1=NULL; // height-filter
	float *mFilterFreArr1=NULL; 
	float *mFilterHeightArr1=NULL;
	int *mIndexArr1=NULL; 
	int *lenArr1=NULL;

	int len=0;
	int len1=0;

	float minHeight=15;

	timeLength=harmonicObj->timeLength;
	peakLength=harmonicObj->peakLength;

	mPeakDbArr=harmonicObj->mPeakDbArr;
	mPeakFreArr=harmonicObj->mPeakFreArr;
	mPeakHeightArr=harmonicObj->mPeakHeightArr;
	mIndexArr=harmonicObj->mIndexArr;
	lenArr=harmonicObj->lenArr;

	mFilterDbArr1=harmonicObj->mFilterDbArr1;
	mFilterFreArr1=harmonicObj->mFilterFreArr1;
	mFilterHeightArr1=harmonicObj->mFilterHeightArr1;
	mIndexArr1=harmonicObj->mIndexArr1;
	lenArr1=harmonicObj->lenArr1;

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

		// fre asc ->start~len1
		__vcorrsort1(mPeakFreArr+(i*peakLength+start),
					mPeakDbArr+(i*peakLength+start), 
					mPeakHeightArr+(i*peakLength+start),
					mIndexArr+(i*peakLength+start),
					len-start, 0);

		for(int j=start;j<len;j++){
			if(mPeakHeightArr[i*peakLength+j]>minHeight){
				float curDb=0;
				float preDb=0;
				float nexDb=0;

				float preHeight=0;
				float nexHeight=0;

				int curIndex=0;
				int preIndex=0;
				int nexIndex=0;

				curDb=mPeakDbArr[i*peakLength+j];
				preDb=mPeakDbArr[i*peakLength+j-1];
				nexDb=mPeakDbArr[i*peakLength+j+1];

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

				if(((curDb-preDb>12)||preHeight>minHeight)&&
					((curDb-nexDb>12)||nexHeight>minHeight)){

					mFilterDbArr1[i*peakLength+len1]=mPeakDbArr[i*peakLength+j];
					mFilterFreArr1[i*peakLength+len1]=mPeakFreArr[i*peakLength+j];
					mFilterHeightArr1[i*peakLength+len1]=mPeakHeightArr[i*peakLength+j];
					mIndexArr1[i*peakLength+len1]=mIndexArr[i*peakLength+j];

					len1++;
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

static void __harmonicObj_filterNear(HarmonicObj harmonicObj){
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

	int lastFlag=1;

	float minFre=30;

	timeLength=harmonicObj->timeLength;
	peakLength=harmonicObj->peakLength;

	mFilterDbArr1=harmonicObj->mFilterDbArr1;
	mFilterFreArr1=harmonicObj->mFilterFreArr1;
	mFilterHeightArr1=harmonicObj->mFilterHeightArr1;
	mIndexArr1=harmonicObj->mIndexArr1;
	lenArr1=harmonicObj->lenArr1;

	mFilterDbArr2=harmonicObj->mFilterDbArr2;
	mFilterFreArr2=harmonicObj->mFilterFreArr2;
	mFilterHeightArr2=harmonicObj->mFilterHeightArr2;
	mIndexArr2=harmonicObj->mIndexArr2;
	lenArr2=harmonicObj->lenArr2;

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

static void __harmonicObj_filterDB(HarmonicObj harmonicObj){
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

	timeLength=harmonicObj->timeLength;
	peakLength=harmonicObj->peakLength;

	mFilterDbArr2=harmonicObj->mFilterDbArr2;
	mFilterFreArr2=harmonicObj->mFilterFreArr2;
	mFilterHeightArr2=harmonicObj->mFilterHeightArr2;
	mIndexArr2=harmonicObj->mIndexArr2;
	lenArr2=harmonicObj->lenArr2;

	mFilterDbArr3=harmonicObj->mFilterDbArr3;
	mFilterFreArr3=harmonicObj->mFilterFreArr3;
	mFilterHeightArr3=harmonicObj->mFilterHeightArr3;
	mIndexArr3=harmonicObj->mIndexArr3;
	lenArr3=harmonicObj->lenArr3;

	maxDBArr=harmonicObj->maxDBArr;

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
				mFilterDbArr3[i*peakLength+j]>-42){

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

void harmonicObj_free(HarmonicObj harmonicObj){

	if(harmonicObj){

		stftObj_free(harmonicObj->stftObj);

		free(harmonicObj->mPowerArr);
		free(harmonicObj->mDbArr);

		free(harmonicObj->mPeakDbArr);
		free(harmonicObj->mPeakFreArr);
		free(harmonicObj->mPeakHeightArr);
		free(harmonicObj->mIndexArr);
		free(harmonicObj->lenArr);

		free(harmonicObj->mFilterDbArr1);
		free(harmonicObj->mFilterFreArr1);
		free(harmonicObj->mFilterHeightArr1);
		free(harmonicObj->mIndexArr1);
		free(harmonicObj->lenArr1);

		free(harmonicObj->mFilterDbArr2);
		free(harmonicObj->mFilterFreArr2);
		free(harmonicObj->mFilterHeightArr2);
		free(harmonicObj->mIndexArr2);
		free(harmonicObj->lenArr2);

		free(harmonicObj->mFilterDbArr3);
		free(harmonicObj->mFilterFreArr3);
		free(harmonicObj->mFilterHeightArr3);
		free(harmonicObj->mIndexArr3);
		free(harmonicObj->lenArr3);

		free(harmonicObj->maxDBArr);

		free(harmonicObj->mRealArr);
		free(harmonicObj->mImageArr);

		free(harmonicObj);
	}
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







