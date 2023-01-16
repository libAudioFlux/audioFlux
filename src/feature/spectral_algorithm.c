// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"

#include "../flux_spectral.h"

#include "spectral_algorithm.h"

struct OpaqueSpectral{
	int num;
	int timeLength;

	float *freBandArr;
	float *energyArr; // timeLength

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

};

// 计算高频sum/c1/c2
static void __spectralObj_calSum(SpectralObj spectralObj,float *mDataArr);
static void __spectralObj_calC1(SpectralObj spectralObj,float *mDataArr);
static void __spectralObj_calC2(SpectralObj spectralObj,float *mDataArr);
static void __spectralObj_calEntropy(SpectralObj spectralObj,float *mDataArr,int isNorm);
static void __spectralObj_calMean(SpectralObj spectralObj,float *mDataArr);

int spectralObj_new(SpectralObj *spectralObj,int num,float *freBandArr){
	int status=0;

	int start=0;
	int end=0;

	int *indexArr=NULL;
	int indexLength=0;

	if(num<2){
		printf("num is error!!!\n");
		return -1;
	}

	SpectralObj sp=NULL;

	sp=*spectralObj=(SpectralObj )calloc(1, sizeof(struct OpaqueSpectral ));

	start=0;
	end=num-1;

	__varangei(start,end+1,1,&indexArr);
	indexLength=num;
	
	sp->num=num;

	sp->freBandArr=freBandArr;

	sp->start=0;
	sp->end=num-1;

	sp->indexArr=indexArr;
	sp->indexLength=indexLength;

	return status;
}

void spectralObj_setTimeLength(SpectralObj spectralObj,int timeLength){
	float *sumArr=NULL;
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;
	float *entropyArr=NULL;

	float *meanFreArr=NULL; 
	float *meanValueArr=NULL; 

	float *energyArr=NULL;

	energyArr=spectralObj->energyArr; // timeLength

	sumArr=spectralObj->sumArr; // timeLength
	cArr1=spectralObj->cArr1;
	cArr2=spectralObj->cArr2;
	entropyArr=spectralObj->entropyArr;

	meanFreArr=spectralObj->meanFreArr;
	meanValueArr=spectralObj->meanValueArr;

	if(spectralObj->timeLength<timeLength||
		spectralObj->timeLength>timeLength*2){
		free(energyArr);

		free(sumArr);
		free(cArr1);
		free(cArr2);
		free(entropyArr);

		free(meanFreArr);
		free(meanValueArr);

		energyArr=__vnew(timeLength, NULL);

		sumArr=__vnew(timeLength, NULL);
		cArr1=__vnew(timeLength, NULL);
		cArr2=__vnew(timeLength, NULL);
		entropyArr=__vnew(timeLength, NULL);

		meanFreArr=__vnew(timeLength, NULL);
		meanValueArr=__vnew(timeLength, NULL);
	}
	
	spectralObj->timeLength=timeLength;

	spectralObj->energyArr=energyArr;

	spectralObj->sumArr=sumArr;
	spectralObj->cArr1=cArr1;
	spectralObj->cArr2=cArr2;
	spectralObj->entropyArr=entropyArr;

	spectralObj->meanFreArr=meanFreArr;
	spectralObj->meanValueArr=meanValueArr;

	// data or memory change, reset
	spectralObj->isSum=0;
	spectralObj->isC1=0;
	spectralObj->isC2=0;
	spectralObj->isEntropy=0;
	spectralObj->isEnNorm=0;

	spectralObj->isMean=0;
}

void spectralObj_setEdge(SpectralObj spectralObj,int start,int end){
	int num=0;

	num=spectralObj->num; // mel/bark/erb num;linear fftLength/2
	if(start>=0&&end<=num-1&&end>start){
		if(start!=spectralObj->start||
			end!=spectralObj->end){

			// 数据变换或内存变换 重置状态
			spectralObj->isSum=0;
			spectralObj->isC1=0;
			spectralObj->isC2=0;
			spectralObj->isEntropy=0;
			spectralObj->isEnNorm=0;

			spectralObj->isMean=0;
		}

		free(spectralObj->indexArr);

		__varangei(start, end+1, 1, &spectralObj->indexArr);
		spectralObj->indexLength=end-start+1;

		spectralObj->start=start;
		spectralObj->end=end;
	}
}

void spectralObj_setEdgeArr(SpectralObj spectralObj,int *indexArr,int indexLength){
	int flag=1;
	int num=0;

	num=spectralObj->num;
	for(int i=0;i<indexLength;i++){
		if(indexArr[i]<0||indexArr[i]>num-1){
			flag=0;
			free(indexArr);
			break;
		}
	}

	if(flag){
		spectralObj->isSum=0;
		spectralObj->isC1=0;
		spectralObj->isC2=0;
		spectralObj->isEntropy=0;
		spectralObj->isEnNorm=0;

		spectralObj->isMean=0;

		free(spectralObj->indexArr);

		spectralObj->indexArr=indexArr;
		spectralObj->indexLength=indexLength;

		spectralObj->start=indexArr[0];
		spectralObj->end=indexArr[indexLength-1];
	}
}

void spectralObj_flatness(SpectralObj spectralObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;
	float *sumArr=NULL;
	
	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;
	sumArr=spectralObj->sumArr;

	if(!spectralObj->isSum){
		__spectralObj_calSum(spectralObj,mDataArr);
	}

	spectral_flatness(mDataArr,nLength,mLength,
					indexArr,indexLength,
					freArr,sumArr,
					dataArr);
}

// isExp 0;step 1,p 1,type 0 sum 1 mean;
void spectralObj_flux(SpectralObj spectralObj,float *mDataArr,
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

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_flux(mDataArr,nLength,mLength,
				indexArr,indexLength,
				step,p,isPostive,_isExp,_type,
				dataArr);
}

void spectralObj_rolloff(SpectralObj spectralObj,float *mDataArr,float threshold,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;
	float *sumArr=NULL;
	
	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;
	sumArr=spectralObj->sumArr;

	if(!spectralObj->isSum){
		__spectralObj_calSum(spectralObj,mDataArr);
	}

	spectral_rolloff(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,sumArr,threshold,
				dataArr);
}

void spectralObj_centroid(SpectralObj spectralObj,float *mDataArr,float *dataArr){

	if(!spectralObj->isC1){
		__spectralObj_calC1(spectralObj,mDataArr);
	}
	
	memcpy(dataArr, spectralObj->cArr1, sizeof(float )*spectralObj->timeLength);
}

void spectralObj_spread(SpectralObj spectralObj,float *mDataArr,float *dataArr){

	if(!spectralObj->isC2){
		__spectralObj_calC2(spectralObj,mDataArr);
	}
	
	memcpy(dataArr, spectralObj->cArr2, sizeof(float )*spectralObj->timeLength);
}

void spectralObj_skewness(SpectralObj spectralObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *sumArr=NULL;
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;
	
	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;

	sumArr=spectralObj->sumArr;
	cArr1=spectralObj->cArr1;
	cArr2=spectralObj->cArr2;

	if(!spectralObj->isC2){
		__spectralObj_calC2(spectralObj,mDataArr);
	}

	spectral_skewness(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,sumArr,cArr1,cArr2,
				dataArr);
}

void spectralObj_kurtosis(SpectralObj spectralObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *sumArr=NULL;
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;
	
	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;

	sumArr=spectralObj->sumArr;
	cArr1=spectralObj->cArr1;
	cArr2=spectralObj->cArr2;

	if(!spectralObj->isC2){
		__spectralObj_calC2(spectralObj,mDataArr);
	}

	spectral_kurtosis(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,sumArr,cArr1,cArr2,
				dataArr);
}

void spectralObj_entropy(SpectralObj spectralObj,float *mDataArr,int isNorm,float *dataArr){

	if(!spectralObj->isEntropy||spectralObj->isEnNorm!=isNorm){
		__spectralObj_calEntropy(spectralObj, mDataArr,isNorm);
	}

	memcpy(dataArr, spectralObj->entropyArr, sizeof(float )*spectralObj->timeLength);
}

void spectralObj_crest(SpectralObj spectralObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;
	float *sumArr=NULL;
	
	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;
	sumArr=spectralObj->sumArr;

	if(!spectralObj->isSum){
		__spectralObj_calSum(spectralObj,mDataArr);
	}

	spectral_crest(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,sumArr,
				dataArr);
}

void spectralObj_slope(SpectralObj spectralObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *meanFreArr=NULL; 
	float *meanValueArr=NULL; 
	
	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;

	meanFreArr=spectralObj->meanFreArr;
	meanValueArr=spectralObj->meanValueArr;

	if(!spectralObj->isMean){
		__spectralObj_calMean(spectralObj,mDataArr);
	}

	spectral_slope(mDataArr,nLength,mLength,
				indexArr,indexLength,
				freArr,meanFreArr,meanValueArr,
				dataArr);
}

void spectralObj_decrease(SpectralObj spectralObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *sumArr=NULL;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	sumArr=spectralObj->sumArr;

	if(!spectralObj->isSum){
		__spectralObj_calSum(spectralObj,mDataArr);
	}

	spectral_decrease(mDataArr,nLength,mLength,
					indexArr,indexLength,
					sumArr,
					dataArr);
}

// p 2 !=0
void spectralObj_bandWidth(SpectralObj spectralObj,float *mDataArr,float p,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;
	float *cArr1=NULL; // μ1
	
	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;
	cArr1=spectralObj->cArr1;

	if(!spectralObj->isC1){
		__spectralObj_calC1(spectralObj,mDataArr);
	}

	spectral_bandWidth(mDataArr,nLength,mLength,
					indexArr,indexLength,
					freArr,cArr1,p,
					dataArr);
}

void spectralObj_rms(SpectralObj spectralObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_rms(mDataArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

// gamma 10
void spectralObj_energy(SpectralObj spectralObj,float *mDataArr,int isLog,float gamma,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	int isPower=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_energy(mDataArr,nLength,mLength,
					indexArr,indexLength,
					isPower,isLog,gamma,
					dataArr);
}

void spectralObj_hfc(SpectralObj spectralObj,float *mDataArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_hfc(mDataArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectralObj_sd(SpectralObj spectralObj,float *mDataArr,int step,int isPostive,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_sd(mDataArr,nLength,mLength,
				indexArr,indexLength,
				step,isPostive,
				dataArr);
}

void spectralObj_sf(SpectralObj spectralObj,float *mDataArr,int step,int isPostive,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_sf(mDataArr,nLength,mLength,
				indexArr,indexLength,
				step,isPostive,
				dataArr);
}

void spectralObj_mkl(SpectralObj spectralObj,float *mDataArr,int type,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_mkl(mDataArr,nLength,mLength,
				indexArr,indexLength,
				type,
				dataArr);
}

void spectralObj_pd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_pd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectralObj_wpd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_wpd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectralObj_nwpd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_nwpd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectralObj_cd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_cd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

void spectralObj_rcd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_rcd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				dataArr);
}

// mag|power threshold >=0
void spectralObj_broadband(SpectralObj spectralObj,float *mDataArr,float threshold,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

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
void spectralObj_novelty(SpectralObj spectralObj,float *mDataArr,
						int step,float threshold,
						SpectralNoveltyMethodType *methodType,SpectralNoveltyDataType *dataType,
						float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	spectral_novelty(mDataArr,nLength,mLength,
					indexArr,indexLength,
					step,threshold,
					methodType,dataType,
					dataArr);
}

void spectralObj_eef(SpectralObj spectralObj,float *mDataArr,int isEnNorm,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *energyArr=NULL; 
	float *entropyArr=NULL;

	float value1=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	energyArr=spectralObj->energyArr;
	entropyArr=spectralObj->entropyArr;

	spectral_energy(mDataArr,nLength,mLength,
					indexArr,indexLength,
					0,0,0,
					energyArr);

	if(!spectralObj->isEntropy||spectralObj->isEnNorm!=isEnNorm){
		__spectralObj_calEntropy(spectralObj,mDataArr,isEnNorm);
	}

	for(int i=0;i<nLength;i++){
		value1=energyArr[i]*entropyArr[i];
		dataArr[i]=sqrtf(1+fabsf(value1));
	}
}

// gamma 1/10/20... song 0.5
void spectralObj_eer(SpectralObj spectralObj,float *mDataArr,int isEnNorm,float gamma,float *dataArr){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *energyArr=NULL; 
	float *entropyArr=NULL;

	float value1=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	energyArr=spectralObj->energyArr;
	entropyArr=spectralObj->entropyArr;

	spectral_energy(mDataArr,nLength,mLength,
					indexArr,indexLength,
					0,0,0,
					energyArr);

	if(!spectralObj->isEntropy||spectralObj->isEnNorm!=isEnNorm){
		__spectralObj_calEntropy(spectralObj,mDataArr,isEnNorm);
	}

	for(int i=0;i<nLength;i++){
		value1=logf(1+energyArr[i]*gamma)/entropyArr[i];
		dataArr[i]=sqrtf(1+fabsf(value1));
	}
}

// statistics
void spectralObj_max(SpectralObj spectralObj,float *mDataArr,float *valueArr,float *freArr2){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float value=0;
	int index=0;
	
	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;
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

void spectralObj_mean(SpectralObj spectralObj,float *mDataArr,float *valueArr,float *freArr2){

	if(!spectralObj->isMean){
		__spectralObj_calMean(spectralObj, mDataArr);
	}

	memcpy(valueArr, spectralObj->meanValueArr, sizeof(float )*spectralObj->timeLength);
	memcpy(freArr2, spectralObj->meanFreArr, sizeof(float )*spectralObj->timeLength);
}

void spectralObj_var(SpectralObj spectralObj,float *mDataArr,float *valueArr,float *freArr2){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *meanValueArr=NULL;
	float *meanFreArr=NULL;

	float value1=0;
	float value2=0;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;

	meanValueArr=spectralObj->meanValueArr;
	meanFreArr=spectralObj->meanFreArr;

	if(indexLength<2){
		return;
	}

	if(!spectralObj->isMean){
		__spectralObj_calMean(spectralObj, mDataArr);
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

// isSum isC1 isC2
static void __spectralObj_calSum(SpectralObj spectralObj,float *mDataArr){
	int num=0;
	int timeLength=0;
	
	int *indexArr=0;
	int indexLength=0;

	float *sumArr=NULL; // timeLength 内存和值都缓存

	num=spectralObj->num;
	timeLength=spectralObj->timeLength;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	sumArr=spectralObj->sumArr;
	for(int i=0;i<timeLength;i++){
		sumArr[i]=0; // 重置数据
		for(int j=0;j<indexLength;j++){
			int _index=0;

			_index=indexArr[j];
			sumArr[i]+=mDataArr[i*num+_index];
		}
	}

	spectralObj->isSum=1;
}

static void __spectralObj_calC1(SpectralObj spectralObj,float *mDataArr){
	int num=0;
	int timeLength=0;
	
	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *sumArr=NULL; // timeLength 内存和值都缓存
	float *cArr1=NULL; // μ1

	num=spectralObj->num;
	timeLength=spectralObj->timeLength;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;

	sumArr=spectralObj->sumArr;
	cArr1=spectralObj->cArr1;

	if(!spectralObj->isSum){
		__spectralObj_calSum(spectralObj,mDataArr);
	}

	spectral_centroid(mDataArr,timeLength,num,
					indexArr,indexLength,
					freArr,sumArr,
					cArr1);

	spectralObj->isSum=1;
	spectralObj->isC1=1;
}

static void __spectralObj_calC2(SpectralObj spectralObj,float *mDataArr){
	int num=0;
	int timeLength=0;
	
	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *sumArr=NULL; // timeLength 内存和值都缓存
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;

	num=spectralObj->num;
	timeLength=spectralObj->timeLength;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;

	sumArr=spectralObj->sumArr;
	cArr1=spectralObj->cArr1;
	cArr2=spectralObj->cArr2;

	if(!spectralObj->isSum){
		__spectralObj_calSum(spectralObj,mDataArr);
	}

	if(!spectralObj->isC1){
		spectral_centroid(mDataArr,timeLength,num,
					indexArr,indexLength,
					freArr,sumArr,
					cArr1);
	}

	spectral_spread(mDataArr,timeLength,num,
					indexArr,indexLength,
					freArr,sumArr,cArr1,
					cArr2);

	spectralObj->isSum=1;
	spectralObj->isC1=1;
	spectralObj->isC2=1;
}

static void __spectralObj_calEntropy(SpectralObj spectralObj,float *mDataArr,int isNorm){
	int nLength=0;
	int mLength=0;

	int *indexArr=0;
	int indexLength=0;

	float *sumArr=NULL;
	float *entropyArr=NULL;

	nLength=spectralObj->timeLength;
	mLength=spectralObj->num;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	sumArr=spectralObj->sumArr;
	entropyArr=spectralObj->entropyArr;

	if(!spectralObj->isSum){
		__spectralObj_calSum(spectralObj,mDataArr);
	}

	spectral_entropy(mDataArr,nLength,mLength,
					indexArr,indexLength,
					sumArr,isNorm,
					entropyArr);

	spectralObj->isEntropy=1;
	spectralObj->isEnNorm=isNorm;
}

static void __spectralObj_calMean(SpectralObj spectralObj,float *mDataArr){
	int num=0;
	int timeLength=0;
	
	int *indexArr=0;
	int indexLength=0;

	float *freArr=NULL;

	float *meanFreArr=NULL; // 仅内存缓存
	float *meanValueArr=NULL; 

	num=spectralObj->num;
	timeLength=spectralObj->timeLength;

	indexArr=spectralObj->indexArr;
	indexLength=spectralObj->indexLength;

	freArr=spectralObj->freBandArr;

	meanFreArr=spectralObj->meanFreArr;
	meanValueArr=spectralObj->meanValueArr;
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

	spectralObj->isMean=1;
}

void spectralObj_free(SpectralObj spectralObj){
	float *freBandArr=NULL;
	float *energyArr=NULL; 

	int *indexArr=NULL;

	float *sumArr=NULL; // timeLength 内存和值都缓存
	float *cArr1=NULL; // μ1
	float *cArr2=NULL;
	float *entropyArr=NULL;

	float *meanFreArr=NULL; // 仅内存缓存
	float *meanValueArr=NULL; 

	if(spectralObj){
		freBandArr=spectralObj->freBandArr;
		energyArr=spectralObj->energyArr;

		indexArr=spectralObj->indexArr;

		sumArr=spectralObj->sumArr;
		cArr1=spectralObj->cArr1;
		cArr2=spectralObj->cArr2;
		entropyArr=spectralObj->entropyArr;

		meanFreArr=spectralObj->meanFreArr;
		meanValueArr=spectralObj->meanValueArr;
	}
}



















