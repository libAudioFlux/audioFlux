// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../dsp/flux_window.h"

#include "auditory_filterBank.h"
#include "nsgt_filterBank.h"

// rect~guass 窗函数设计法
static void __nsgt_standardFilterBank(int num,
									SpectralFilterBankStyleType styleType,
									SpectralFilterBankNormalType normType,
									int *lenArr,int *binBandArr,
									float *windowDataArr,int *indexArr);

static void __nsgt_efficientFilterBank1(int num,
									SpectralFilterBankStyleType styleType,
									SpectralFilterBankNormalType normType,
									int *lenArr,int *binBandArr,
									float *windowDataArr,int *indexArr);

static void __nsgt_efficientFilterBank2(int num,
									SpectralFilterBankStyleType styleType,
									SpectralFilterBankNormalType normType,
									int *lenArr,int *binBandArr,
									float *windowDataArr,int *indexArr);

static void __nsgt_calBandEdge(int num,int dataLength,int samplate,
							float lowFre,float highFre,int isEdge,
							void *callback1,void *callback2,float ref,
							float **freBandArr,int **binBandArr);

/***
	nsgt-filterBank相关 
		log+linear+mel/bark/erb 都不包含边界
			log 针对nsgt-spectrogram形式上包含
		log is nsg-cqt
	SpectralFilterBankScaleType 默认linear linear~Log
	SpectralFilterBankStyleType 默认Slaney 除Gammatone
	SpectralFilterBankNormalType 默认None 除area
****/
void nsgt_filterBank(int num,int dataLength,int samplate, int minLength,
					int type,
					SpectralFilterBankScaleType scaleType,
					SpectralFilterBankStyleType styleType,
					SpectralFilterBankNormalType normType,
					float lowFre,float highFre,int binPerOctave,
					float **filterBankArr,int *lengthArr,
					float *freBandArr,int *binBandArr,int *offsetArr,
					int *maxLength,int *totalLength){
	int totalLen=0;

	float *fArr=NULL;
	int *bArr=NULL;
	int *iArr=NULL;

	int isEdge=0;
	int offset=0;

	float *winDataArr=NULL;
	int *lenArr=NULL;

	// UniFunc||UniFunc1
	void *func1=NULL;
	void *func2=NULL;

	float ref=0;

	if(styleType==SpectralFilterBankStyle_Gammatone){ // include edge!!!
		isEdge=1;
	}

	if(!isEdge){
		offset=1;
	}

	// 0. revise log/linear!!!
	if(scaleType==SpectralFilterBankScale_Octave){ 
		if(binPerOctave&&binPerOctave>=4&&binPerOctave<=48){
			ref=binPerOctave;
		}
		else{
			ref=12;
		}

		auditory_reviseLogFre(num,lowFre,highFre,ref,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Linear){
		ref=samplate*1.0/dataLength;

		auditory_reviseLinearFre(num,lowFre,highFre,ref,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Linspace){
		auditory_reviseLinspaceFre(num,lowFre,highFre,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Log){
		auditory_reviseLogspaceFre(num,lowFre,highFre,isEdge,&lowFre,&highFre);
	}

	// 1. cal freBand/binBand
	if(scaleType==SpectralFilterBankScale_Linear){ // linear
		func1=auditory_freToLinear;
		func2=auditory_linearToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Linspace){ // linspace
		func1=auditory_freToLinspace;
		func2=auditory_linspaceToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Mel){ // mel
		func1=auditory_freToMel;
		func2=auditory_melToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Bark){ // bark
		func1=auditory_freToBark;
		func2=auditory_barkToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Erb){ // erb
		func1=auditory_freToErb;
		func2=auditory_erbToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Octave){ // log
		func1=auditory_freToLog;
		func2=auditory_logToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Log){ // logspace
		func1=auditory_freToLogspace;
		func2=auditory_logspaceToFre;
	}
	
	__nsgt_calBandEdge(num,dataLength,samplate,
					lowFre,highFre,isEdge,
					func1,func2,ref,
					&fArr,&bArr);

	// 2. 计算filterBank
	lenArr=__vnewi(num, NULL);
	__vsubi(bArr+2, bArr, num, lenArr);

	if(type){ // 1 NSGTFilterBank_Standard "periodic" even or odd
		for(int i=0;i<num;i++){ // width+1 == order
			lenArr[i]+=1;
			if(lenArr[i]<minLength){
				lenArr[i]=minLength;
			}
		}
	}
	else{ // 0 NSGTFilterBank_Efficient "symmetric"
		for(int i=0;i<num;i++){ // width+1 == order
			int left=0;
			int cur=0;
			int right=0;

			int v1=0;
			int v2=0;
			int v3=0;

			left=bArr[i];
			cur=bArr[i+1];
			right=bArr[i+2];

			if(right-left>=1){ // >= 3 point
				v1=cur-left;
				v2=right-cur;

				v3=(v2>=v1?v2:v1);
				lenArr[i]=2*v3+1;
			}
			else{
				lenArr[i]=0;
			}
			
			if(lenArr[i]<minLength){
				lenArr[i]=minLength;
			}
		}
	}

	totalLen=__vsumi(lenArr, num);
	winDataArr=__vnew(2*totalLen, NULL);
	iArr=__vnewi(num, NULL);

	if(!type){ // efficient 1 max 2 suitable
		__nsgt_efficientFilterBank1(num,
								styleType,
								normType,
								lenArr,bArr,
								winDataArr,iArr);
	}
	else{ // standard
		__nsgt_standardFilterBank(num,
								styleType,
								normType,
								lenArr,bArr,
								winDataArr,iArr);
	}

	// 3. 处理freBand/binBand
	if(freBandArr){
		memcpy(freBandArr, fArr+offset, sizeof(float )*num);
	}
	
	if(binBandArr){
		memcpy(binBandArr, bArr+offset, sizeof(int )*num);
	}

	if(offsetArr){
		memcpy(offsetArr, iArr, sizeof(int )*num);
	}

	if(lengthArr){
		memcpy(lengthArr, lenArr, sizeof(int )*num);
	}

	if(filterBankArr){
		*filterBankArr=winDataArr;
	}
	else{
		free(winDataArr);
	}

	if(maxLength){
		__vmaxi(lenArr, num, maxLength);
	}

	if(totalLength){
		*totalLength=totalLen;
	}

	free(fArr);
	free(bArr);
	free(lenArr);
	free(iArr);
}

/***
	lenArr 一般有最小限制 minWindow>=2
	标准论文直接取len长度的window数据,以当前bin为中心参与后续计算
	实际上,当前bin前后带宽通常情况下并不一致,一般后面大于前面
	导致每个频带都会有之前累加,低频区域会出现明显叠影
****/
static void __nsgt_standardFilterBank(int num,
									SpectralFilterBankStyleType styleType,
									SpectralFilterBankNormalType normType,
									int *lenArr,int *binBandArr,
									float *windowDataArr,int *offsetArr){
	int len=0;
	int index=0;

	for(int i=0;i<num;i++){
		float *arr=NULL;

		len=lenArr[i];
		offsetArr[i]=binBandArr[i+1]-len/2;

		if(offsetArr[i]<0){
			offsetArr[i]=0;
		}

		if(styleType==SpectralFilterBankStyle_Slaney){ // Triang
			arr=window_createTriang(len, 1);
		}
		else if(styleType==SpectralFilterBankStyle_ETSI){ // Bartlett
			arr=window_createBartlett(len, 1);
		}
		else if(styleType==SpectralFilterBankStyle_Hann){
			arr=window_createHann(len, 1);
		}
		else if(styleType==SpectralFilterBankStyle_Hamm){
			arr=window_createHamm(len, 1);
		}
		else if(styleType==SpectralFilterBankStyle_Blackman){
			arr=window_createBlackman(len, 1);
		}
		else if(styleType==SpectralFilterBankStyle_Bohman){
			arr=window_createBohman(len, 1);
		}
		else if(styleType==SpectralFilterBankStyle_Kaiser){
			arr=window_createKaiser(len, 1,NULL);
		}
		else if(styleType==SpectralFilterBankStyle_Gauss){
			arr=window_createGauss(len, 1,NULL);
		}
		else{ // Rect
			arr=__vnew(len, NULL);
			for(int i=0;i<len;i++){
				arr[i]=1;
			}
		}

		if(normType==SpectralFilterBankNormal_BandWidth){
			__vdiv_value(arr, sqrtf(len), len, NULL); // sqrt(len) or len/2 ???
		}

		memcpy(windowDataArr+index, arr, sizeof(float )*len);
		index+=len;

		free(arr);
	}
}

static void __nsgt_efficientFilterBank1(int num,
									SpectralFilterBankStyleType styleType,
									SpectralFilterBankNormalType normType,
									int *lenArr,int *binBandArr,
									float *windowDataArr,int *offsetArr){
	int len=0;
	int index=0;

	for(int i=0;i<num;i++){
		float *arr=NULL;

		len=lenArr[i];
		offsetArr[i]=binBandArr[i+1]-len/2;

		if(offsetArr[i]<0){
			offsetArr[i]=0;
		}

		if(styleType==SpectralFilterBankStyle_Slaney){ // Triang
			arr=window_createTriang(len, 0);
		}
		else if(styleType==SpectralFilterBankStyle_ETSI){ // Bartlett
			arr=window_createBartlett(len, 0);
		}
		else if(styleType==SpectralFilterBankStyle_Hann){
			arr=window_createHann(len, 0);
		}
		else if(styleType==SpectralFilterBankStyle_Hamm){
			arr=window_createHamm(len, 0);
		}
		else if(styleType==SpectralFilterBankStyle_Blackman){
			arr=window_createBlackman(len, 0);
		}
		else if(styleType==SpectralFilterBankStyle_Bohman){
			arr=window_createBohman(len, 0);
		}
		else if(styleType==SpectralFilterBankStyle_Kaiser){
			arr=window_createKaiser(len, 0,NULL);
		}
		else if(styleType==SpectralFilterBankStyle_Gauss){
			arr=window_createGauss(len, 0,NULL);
		}
		else{ // Rect
			arr=__vnew(len, NULL);
			for(int i=0;i<len;i++){
				arr[i]=1;
			}
		}

		if(normType==SpectralFilterBankNormal_BandWidth){
			__vdiv_value(arr, sqrtf(len), len, NULL); // sqrt(len) or len/2 ???
		}

		memcpy(windowDataArr+index, arr, sizeof(float )*len);
		index+=len;

		free(arr);
	}
}

static void __nsgt_efficientFilterBank2(int num,
									SpectralFilterBankStyleType styleType,
									SpectralFilterBankNormalType normType,
									int *lenArr,int *binBandArr,
									float *windowDataArr,int *offsetArr){
	int len=0;
	int index=0;

	int left=0;
	int cur=0;
	int right=0;

	for(int i=0;i<num;i++){
		len=lenArr[i];

		left=binBandArr[i];
		cur=binBandArr[i+1];
		right=binBandArr[i+2];

		if(cur>left){
			float *wArr=NULL;
			int det=0;
			int offset=0;

			det=cur-left;
			offsetArr[i]=left;
			
			// if(offsetArr[i]<0){
			// 	offsetArr[i]=0;
			// }

			if(styleType==SpectralFilterBankStyle_Hann){
				wArr=window_createHann(2*det+1,0);
			}
			else if(styleType==SpectralFilterBankStyle_Hamm){
				wArr=window_createHamm(2*det+1,0);
			}
			else if(styleType==SpectralFilterBankStyle_Blackman){
				wArr=window_createBlackman(2*det+1,0);
			}
			else if(styleType==SpectralFilterBankStyle_Bohman){
				wArr=window_createBohman(2*det+1,0);
			}
			else if(styleType==SpectralFilterBankStyle_Kaiser){
				wArr=window_createKaiser(2*det+1,0,NULL);
			}
			else if(styleType==SpectralFilterBankStyle_Gauss){ // gauss
				wArr=window_createGauss(2*det+1,0,NULL);
			}
			else{ // Rect
				wArr=__vnew(2*det+1, NULL);
				for(int i=0;i<2*det+1;i++){
					wArr[i]=1;
				}
			}

			if(normType==SpectralFilterBankNormal_BandWidth){
				__vdiv_value(wArr, sqrtf(len), 2*det+1, NULL); // sqrt(len) or len/2 ???
			}

			memcpy(windowDataArr+index, wArr+offset, sizeof(float )*(det+1));
			index+=(det+1);

			free(wArr);
		}

		if(right>cur){
			float *wArr=NULL;
			int det=0;
			int offset=0;

			det=right-cur;
			offset=(2*det+1)/2+1;
			if(styleType==SpectralFilterBankStyle_Hann){
				wArr=window_createHann(2*det+1,0);
			}
			else if(styleType==SpectralFilterBankStyle_Hamm){
				wArr=window_createHamm(2*det+1,0);
			}
			else if(styleType==SpectralFilterBankStyle_Blackman){
				wArr=window_createBlackman(2*det+1,0);
			}
			else if(styleType==SpectralFilterBankStyle_Bohman){
				wArr=window_createBohman(2*det+1,0);
			}
			else if(styleType==SpectralFilterBankStyle_Kaiser){
				wArr=window_createKaiser(2*det+1,0,NULL);
			}
			else if(styleType==SpectralFilterBankStyle_Gauss){ // gauss
				wArr=window_createGauss(2*det+1,0,NULL);
			}
			else{ // Rect
				wArr=__vnew(2*det+1, NULL);
				for(int i=0;i<2*det+1;i++){
					wArr[i]=1;
				}
			}

			if(normType==SpectralFilterBankNormal_BandWidth){
				__vdiv_value(wArr, sqrtf(len), 2*det+1, NULL); // sqrt(len) or len/2 ???
			}

			memcpy(windowDataArr+index, wArr+offset, sizeof(float )*det);
			index+=det;

			free(wArr);
		}

		if(left==right){ // only 1 set 0||1 ???
			windowDataArr[index]=1; // 0 ???
			index++;
		}
	}
}

static void __nsgt_calBandEdge(int num,int dataLength,int samplate,
							float lowFre,float highFre,int isEdge,
							void *callback1,void *callback2,float ref,
							float **freBandArr,int **binBandArr){
	float *fArr=NULL;
	int *bArr=NULL;

	float low=0;
	float high=0;

	int det=0;

	UniFunc func1=NULL;
	UniFunc func2=NULL;

	UniFunc1 func3=NULL;
	UniFunc1 func4=NULL;

	if(!isEdge){ // ! gammatone&linear
		det=2;
	}

	if(callback1&&callback2){
		if(!ref){
			func1=(UniFunc )callback1;
			func2=(UniFunc )callback2;

			low=func1(lowFre);
			high=func1(highFre);
		}
		else{
			func3=(UniFunc1 )callback1;
			func4=(UniFunc1 )callback2;

			low=func3(lowFre,ref);
			high=func3(highFre,ref);
		}
	}
	else{
		low=lowFre;
		high=highFre;
	}

	// 1. mel/barkArr
	fArr=__vlinspace(low, high, num+det, 0);

	// 2. mel/barkArr=>freArr
	if(callback1&&callback2){
		if(!ref){
			__vmap(fArr, num+det, func2, NULL);
		}
		else{
			__vmap1(fArr, num+det, func4,ref, NULL);
		}
	}
	
	if(freBandArr){
		*freBandArr=fArr;
	}

	// 3. freArr=>binArr vector???
	if(binBandArr){
		bArr=__vnewi(num+det, NULL);
		for(int i=0;i<num+det;i++){
			bArr[i]=roundf(dataLength*fArr[i]/samplate);
		}

		*binBandArr=bArr;
	}

	if(!freBandArr){
		free(fArr);
	}
}










