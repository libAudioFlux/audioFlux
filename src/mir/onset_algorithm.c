// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "../util/flux_util.h"

#include "../flux_spectral.h"

#include "onset_algorithm.h"

struct OpaqueOnset{
	NoveltyType noveltyType;

	int nLength; // timeLength
	int mLength; // freLength

	int order; // fre

	int preMax;
	int postMax;

	int preAvg;
	int postAvg;
	
	int wait;
	float delta;

	float *mFilterArr; // nLength*mLength
	float *mFluxArr; // (nLength-step)*mLength

	// novelty params
	int step;
	float p;
	int isPostive;
	int isExp;
	int type;

	float threshold;

	int isNorm;
	float gamma;
	
};

static void _onsetObj_dealFilterArr(OnsetObj onsetObj,float *mDataArr1);
static void _onsetObj_dealFluxArr(OnsetObj onsetObj,float *mDataArr1,float *mDataArr2,int *indexArr,int indexLength,float *vArr1);
static int _onsetObj_dealPointArr(OnsetObj onsetObj,float *evnArr,int *pointArr);

static void _onsetObj_initPeak(OnsetObj onsetObj,int samplate,int slideLength);
static void _onsetObj_initParam(OnsetObj onsetObj,NoveltyParam *param);

static int __peakPick(float *evnArr,int *pointArr,int length,int preMax,int postMax,int preAvg,int postAvg,int wait,float delta);

int onsetObj_new(OnsetObj *onsetObj,int nLength,int mLength,int slideLength,
				int *samplate,int *filterOrder,
				NoveltyType *type){
	int status=0;

	NoveltyType _type=Novelty_Flux;

	int order=1;

	float *mFluxArr=NULL;
	float *mFilterArr=NULL;

	int _samplate=32000;

	OnsetObj onset=NULL;

	onset=*onsetObj=(OnsetObj )calloc(1,sizeof(struct OpaqueOnset ));

	if(type){
		_type=*type;
	}

	if(filterOrder){
		if(*filterOrder>0){
			order=*filterOrder;
		}
	}

	if(samplate){
		if(*samplate>0){
			_samplate=*samplate;
		}
	}

	if(slideLength<1){
		slideLength=512;
	}

	mFluxArr=__vnew(nLength*mLength, NULL);
	mFilterArr=__vnew(nLength*mLength, NULL);

	onset->noveltyType=_type;

	onset->nLength=nLength;
	onset->mLength=mLength;

	onset->order=order;

	onset->mFluxArr=mFluxArr;
	onset->mFilterArr=mFilterArr;
	
	_onsetObj_initPeak(onset,_samplate,slideLength);

	return status;
}

/***
	large-scale from CPJKU onset_db
	preMax 0.03*sr/slideLength
	postMax 0.00*sr/slideLength+1
	preAvg 0.1*sr/slideLength
	postAvg 0.1*sr/slideLength+1
	wait 0.03*sr/slideLength
	delta 0.07
****/
static void _onsetObj_initPeak(OnsetObj onsetObj,int samplate,int slideLength){

	onsetObj->preMax=floorf(0.03*samplate/slideLength);
	onsetObj->postMax=floorf(0.0*samplate/slideLength+1);

	onsetObj->preAvg=floorf(0.1*samplate/slideLength);
	onsetObj->postAvg=floorf(0.1*samplate/slideLength+1);

	onsetObj->wait=floorf(0.03*samplate/slideLength);
	onsetObj->delta=0.07;
}

static void _onsetObj_initParam(OnsetObj onsetObj,NoveltyParam *param){
	int step=1;
	float p=1;
	int isPostive=1;
	int isExp=0;
	int type=0;

	float threshold=0;

	int isNorm=0;
	float gamma=1;

	if(param){
		if(param->step>0){
			step=param->step;
		}

		if(param->p){ // !=0
			p=param->p;
		}

		isPostive=param->isPostive;
		isExp=param->isExp;
		type=param->type;

		threshold=param->threshold;

		isNorm=param->isNorm;
		if(param->gamma>0){
			gamma=param->gamma;
		}
	}


	onsetObj->step=step;
	onsetObj->p=p;
	onsetObj->isPostive=isPostive;
	onsetObj->isExp=isExp;
	onsetObj->type=type;

	onsetObj->threshold=threshold;

	onsetObj->isNorm=isNorm;
	onsetObj->gamma=gamma;
}

/***
	1. cal evnArr
	2. cal pointArr
****/
int onsetObj_onset(OnsetObj onsetObj,float *mDataArr1,float *mDataArr2,
				NoveltyParam *param,int *indexArr,int indexLength,
				float *evnArr,int *pointArr){
	int pointLength=0;

	int *iArr=NULL;
	int iLength=0;

	if(indexArr){
		iArr=__vnewi(indexLength, NULL);
		memcpy(iArr, indexArr, sizeof(int )*indexLength);
		iLength=indexLength;
	}
	else{
		__varangei(0, onsetObj->mLength, 1, &iArr);
		iLength=onsetObj->mLength;
	}

	_onsetObj_initParam(onsetObj, param);

	_onsetObj_dealFilterArr(onsetObj,mDataArr1);
	_onsetObj_dealFluxArr(onsetObj,mDataArr1,mDataArr2,iArr,iLength,evnArr);
	pointLength=_onsetObj_dealPointArr(onsetObj,evnArr,pointArr);

	free(iArr);
	return pointLength;
}

static int _onsetObj_dealPointArr(OnsetObj onsetObj,float *evnArr,int *pointArr){
	int pointLength=0;
	int length=0;

	int preMax=0;
	int postMax=0;
	int preAvg=0;
	int postAvg=0;
	int wait=0;
	float delta=0;

	length=onsetObj->nLength;

	preMax=onsetObj->preMax;
	postMax=onsetObj->postMax;
	preAvg=onsetObj->preAvg;
	postAvg=onsetObj->postAvg;
	delta=onsetObj->delta;
	wait=onsetObj->wait;

	pointLength=__peakPick(evnArr, pointArr, length, preMax, postMax, preAvg, postAvg,wait,delta);

	return pointLength;
}

static void _onsetObj_dealFilterArr(OnsetObj onsetObj,float *mDataArr1){
	int nLength=0; // timeLength
	int mLength=0; // freLength

	int order=0; // fre

	float *mFilterArr=NULL; // nLength*mLength

	order=onsetObj->order;
	if(order<2){
		return;
	}

	nLength=onsetObj->nLength;
	mLength=onsetObj->mLength;

	mFilterArr=onsetObj->mFilterArr;

	__mmaxfilter(mDataArr1,nLength,mLength,1,order,mFilterArr);
}

static void _onsetObj_dealFluxArr(OnsetObj onsetObj,float *mDataArr1,float *mDataArr2,int *indexArr,int indexLength,float *vArr1){
	NoveltyType noveltyType=Novelty_Flux;

	int nLength=0; // timeLength
	int mLength=0; // freLength

	int order=0; // fre

	float *mFluxArr=NULL; // (nLength-step)*mLength
	float *mArr=NULL;

	float min=0;
	float max=0;

	int step=1; // >=1
	float p=1; // >=1
	int isPostive=1; // 1
	int isExp=0; // 0
	int type=0; // 0 sum 1 mean

	float threshold=0; // >=0

	int isNorm=0; // 0|1
	float gamma=1;

	int start=0;
	int end=0;

	noveltyType=onsetObj->noveltyType;

	nLength=onsetObj->nLength;
	mLength=onsetObj->mLength;

	order=onsetObj->order;

	mFluxArr=onsetObj->mFluxArr;
	if(order>1){
		mArr=onsetObj->mFilterArr;
	}
	else{
		mArr=mDataArr1;
	}

	start=0;
	end=nLength-1;

	step=onsetObj->step;
	p=onsetObj->p;
	isPostive=onsetObj->isPostive;
	isExp=onsetObj->isExp;
	type=onsetObj->type;

	threshold=onsetObj->threshold;

	isNorm=onsetObj->isNorm;
	gamma=onsetObj->gamma;

	memset(vArr1, 0, sizeof(float )*step);
	
	if(noveltyType==Novelty_HFC){
		spectral_hfc(mArr,nLength,mLength,
					indexArr,indexLength,
					vArr1);
	}
	else if(noveltyType==Novelty_SD){
		spectral_sd(mArr,nLength,mLength,
					indexArr,indexLength,
					step,isPostive,
					vArr1);
	}
	else if(noveltyType==Novelty_SF){
		spectral_sf(mArr,nLength,mLength,
					indexArr,indexLength,
					step,isPostive,
					vArr1);
	}
	else if(noveltyType==Novelty_MKL){
		spectral_mkl(mArr,nLength,mLength,
					indexArr,indexLength,
					type,
					vArr1);
	}
	else if(noveltyType==Novelty_PD){
		spectral_pd(mArr,mDataArr2,nLength,mLength,
					indexArr,indexLength,
					vArr1);
	}
	else if(noveltyType==Novelty_WPD){
		spectral_wpd(mArr,mDataArr2,nLength,mLength,
					indexArr,indexLength,
					vArr1);
	}
	else if(noveltyType==Novelty_NWPD){
		spectral_nwpd(mArr,mDataArr2,nLength,mLength,
					indexArr,indexLength,
					vArr1);
	}
	else if(noveltyType==Novelty_CD){
		spectral_cd(mArr,mDataArr2,nLength,mLength,
					indexArr,indexLength,
					vArr1);
	}
	else if(noveltyType==Novelty_RCD){
		spectral_rcd(mArr,mDataArr2,nLength,mLength,
					indexArr,indexLength,
					vArr1);
	}
	else if(noveltyType==Novelty_Broadband){
		spectral_broadband(mArr,nLength,mLength,
							indexArr,indexLength,
							threshold,
							vArr1);
	}
	else{ // flux
		spectral_flux(mArr,nLength,mLength,
					indexArr,indexLength,
					step,p,isPostive,isExp,type,
					vArr1);
	}

	__vmin(vArr1, nLength, &min);
	__vsub_value(vArr1, min, nLength, vArr1);

	__vmax(vArr1, nLength, &max);
	if(max>0){
		__vdiv_value(vArr1, max, nLength, vArr1);
	}
}

void onsetObj_free(OnsetObj onsetObj){
	float *mFilterArr=NULL; // nLength*mLength
	float *mFluxArr=NULL; // (nLength-step)*mLength

	if(!onsetObj){
		return;
	}

	mFilterArr=onsetObj->mFilterArr;
	mFluxArr=onsetObj->mFluxArr;

	free(mFilterArr);
	free(mFluxArr);

	free(onsetObj);
}

void onsetObj_debug(OnsetObj onsetObj){

	printf("onsetObj is :\n");
	printf("preMax=%d,postMax=%d, preAvg=%d,postAvg=%d, wait=%d,delta=%f\n",
			onsetObj->preMax,onsetObj->postMax,onsetObj->preAvg,onsetObj->postAvg,onsetObj->wait,onsetObj->delta);

	printf("timeLength=%d,freNum=%d, step=%d,order=%d\n",
			onsetObj->nLength,onsetObj->mLength,onsetObj->step,onsetObj->order);

	printf("\n");
}

// peak
/***
	x[n]==max(x[n-preMax:n+postMax]) &&
	x[n]>=mean(x[n-preAvg:n+postAvg])+delta &&
	n-pre>wait => n is point
****/
static int __peakPick(float *evnArr,int *pointArr,int length,int preMax,int postMax,int preAvg,int postAvg,int wait,float delta){
	int pointLength=0;

	int pre=0;

	float max=0;
	float mean=0;

	int start1=0;
	int end1=0;

	int start2=0;
	int end2=0;

	pre=-wait-1;
	for(int i=0;i<length;i++){
		start1=(i-preMax>=0?i-preMax:0);
		end1=(i+postMax<length?i-1+postMax:length-1);

		__vmax(evnArr+start1, end1-start1+1, &max);
		if(evnArr[i]==max){
			start2=(i-preAvg>=0?i-preAvg:0);
			end2=(i+postAvg<length?i-1+postAvg:length-1);

			mean=__vmean(evnArr+start2, end2-start2+1);
			if(evnArr[i]>=mean+delta){
				if(i-pre>wait){
					pointArr[pointLength]=i;
					pre=i;

					pointLength++;
				}
			}
		}
	}

	return pointLength;
}

