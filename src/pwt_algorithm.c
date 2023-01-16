// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/flux_correct.h"
#include "dsp/fft_algorithm.h"
#include "dsp/dft_algorithm.h"

#include "filterbank/auditory_filterBank.h"

#include "pwt_algorithm.h"

struct OpaquePWT{
	FFTObj fftObj;
	DFTObj dftObj;

	int fftLength; // 1<<radix2Exp+2*padLength
	int dataLength; // data length=1<<radix2Exp
	int padLength; // 0||dataLength/2||log2(dataLength) ???

	int num;

	int *binBandArr;
	float *freBandArr;

	float *mFilterBankArr; // num*fftLength
	float *mFilterBankDetArr; // j

	float *realArr1; // fft data result
	float *imageArr1;

	float *mRealArr2; // mFilterBankArr.dot(fftData)  
	float *mImageArr2;

	float *mRealArr3; // ifft(mArr2) !isPad not use
	float *mImageArr3;

	float *curDataArr; // cache reflect padd

	// params
	int samplate; // filterBank 
	float lowFre;
	float highFre;

	int binPerOctave; 

	SpectralFilterBankScaleType scaleType;
	SpectralFilterBankStyleType styleType;
	SpectralFilterBankNormalType normalType;

};

static void __pwtObj_init(PWTObj pwtObj);
static void __pwtObj_pwt(PWTObj pwtObj,float *dataArr,float *mFilterBankArr,int iFlag,float *mRealArr4,float *mImageArr4);

int pwtObj_new(PWTObj *pwtObj,int num,int radix2Exp,
			 int *samplate,float *lowFre,float *highFre,int *binPerOctave,
			 SpectralFilterBankScaleType *scaleType,
			 SpectralFilterBankStyleType *styleType,
			 SpectralFilterBankNormalType *normalType,
			 int *isPadding){
	int status=0;
	PWTObj pwt=NULL;

	FFTObj fftObj=NULL;
	DFTObj dftObj=NULL;

	int fftLength=0; // 1<<radix2Exp+2*padLength
	int dataLength=0; // data length=1<<radix2Exp
	int padLength=0;

	float *realArr1=NULL; // fft data result
	float *imageArr1=NULL;

	float *mRealArr2=NULL; // mFilterBankArr.dot(fftData)  
	float *mImageArr2=NULL;

	float *mRealArr3=NULL; // ifft(mArr2) !isPad not use
	float *mImageArr3=NULL;

	float *curDataArr=NULL;

	int _samplate=32000; // filterBank 
	float _lowFre=0;
	float _highFre=0;

	int _binPerOctave=12; // log

	int _isPadding=0;

	int _flag=0;

	SpectralFilterBankScaleType _scaleType=SpectralFilterBankScale_Octave;
	SpectralFilterBankStyleType _styleType=SpectralFilterBankStyle_Slaney;
	SpectralFilterBankNormalType _normalType=SpectralFilterBankNormal_None;

	if(radix2Exp){
		if(radix2Exp<1||radix2Exp>30){
			status=-100;
			printf("radix2Exp is error!\n");
			return status;
		}
	}

	fftLength=1<<radix2Exp;
	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
		}
	}

	if(scaleType){
		_scaleType=*scaleType;
		if(_scaleType>SpectralFilterBankScale_Log){
			printf("scaleType is error!\n");
			return 1;
		}
	}

	if(styleType){
		_styleType=*styleType;
	}

	if(normalType){
		_normalType=*normalType;
	}

	_highFre=_samplate/2.0;

	if(lowFre){
		if(*lowFre>=0&&*lowFre<_samplate/2.0){
			_lowFre=*lowFre;
		}
	}

	if(_lowFre==0){
		if(_scaleType==SpectralFilterBankScale_Octave||
			_scaleType==SpectralFilterBankScale_Log){ // Log/Logspace

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
		if(_scaleType==SpectralFilterBankScale_Octave||
			_scaleType==SpectralFilterBankScale_Log){ // Log/Logspace
			
			_lowFre=powf(2, -45/12.0)*440;
			_highFre=powf(2, 38/12.0)*440;
		}
	}

	if(binPerOctave){
		if(*binPerOctave>=4&&*binPerOctave<=48){
			_binPerOctave=*binPerOctave;
		}
	}

	if(_scaleType==SpectralFilterBankScale_Linear){ // linear
		float det=0;

		det=_samplate/(float )fftLength;
		auditory_reviseLinearFre(num, _lowFre,_highFre, det,1,&_lowFre, &_highFre);

		if(_highFre>_samplate/2.0){
			printf("scale linear: lowFre and num is large, overflow error\n");
			return -1;
		}
	}
	else if(_scaleType==SpectralFilterBankScale_Octave){ // log
		auditory_reviseLogFre(num, _lowFre,_highFre, _binPerOctave,1, &_lowFre, &_highFre);

		if(_highFre>_samplate/2.0){
			printf("scale log: lowFre and num is large, overflow error!\n");
			return -1;
		}
	}

	if(isPadding){
		_isPadding=*isPadding;
	}

	if(num<2||num>fftLength/2+1){ 
		printf("num is error!\n");
		return -1;
	}

	pwt=*pwtObj=(PWTObj )calloc(1, sizeof(struct OpaquePWT ));

	dataLength=1<<radix2Exp;
	if(_isPadding){
		if(dataLength<=1e5){
			padLength=dataLength/2;
		}
		else{
			padLength=ceilf(log2f(dataLength));
		}

		fftLength=dataLength+2*padLength;

		curDataArr=__vnew(fftLength, NULL);
		mRealArr3=__vnew(num*fftLength, NULL);
		mImageArr3=__vnew(num*fftLength, NULL);
	}
	else{
		fftLength=dataLength;
	}

	_flag=util_isPowerTwo(fftLength);
	if(_flag){
		radix2Exp=util_powerTwoBit(fftLength);
		fftObj_new(&fftObj, radix2Exp);
	}
	else{
		dftObj_new(&dftObj, fftLength);
	}

	realArr1=__vnew(fftLength, NULL);
	imageArr1=__vnew(fftLength, NULL);

	mRealArr2=__vnew(num*fftLength, NULL);
	mImageArr2=__vnew(num*fftLength, NULL);

	pwt->fftObj=fftObj;
	pwt->dftObj=dftObj;

	pwt->fftLength=fftLength;
	pwt->dataLength=dataLength;
	pwt->padLength=padLength;

	pwt->num=num;

	pwt->realArr1=realArr1;
	pwt->imageArr1=imageArr1;

	pwt->mRealArr2=mRealArr2;
	pwt->mImageArr2=mImageArr2;

	pwt->mRealArr3=mRealArr3;
	pwt->mImageArr3=mImageArr3;

	pwt->fftObj=fftObj;

	pwt->curDataArr=curDataArr;

	pwt->samplate=_samplate;
	pwt->lowFre=_lowFre;
	pwt->highFre=_highFre;

	pwt->binPerOctave=_binPerOctave;

	pwt->scaleType=_scaleType;
	pwt->styleType=_styleType;
	pwt->normalType=_normalType;

	__pwtObj_init(pwt);

	return status;
}

static void __pwtObj_init(PWTObj pwtObj){
	int *binBandArr=NULL;
	float *freBandArr=NULL;

	float *mFilterBankArr=NULL;

	int num=0;
	int dataLength=0;
	int samplate=0;
	int padLength=0;
	int fftLength=0;

	SpectralFilterBankScaleType scaleType;
	SpectralFilterBankStyleType styleType;
	SpectralFilterBankNormalType normalType;

	float lowFre=0;
	float highFre=0;
	int binPerOctave=0;

	num=pwtObj->num;
	dataLength=pwtObj->dataLength;
	samplate=pwtObj->samplate;
	padLength=pwtObj->padLength;
	fftLength=pwtObj->fftLength;

	scaleType=pwtObj->scaleType;
	styleType=pwtObj->styleType;
	normalType=pwtObj->normalType;

	lowFre=pwtObj->lowFre;
	highFre=pwtObj->highFre;
	binPerOctave=pwtObj->binPerOctave;

	mFilterBankArr=__vnew(num*fftLength, NULL);
	freBandArr=__vnew(num+2, NULL);
	binBandArr=__vnewi(num+2, NULL);

	auditory_filterBank(num,fftLength,samplate,1,
						scaleType,styleType,normalType,
						lowFre,highFre,binPerOctave,
						mFilterBankArr,
						freBandArr,
						binBandArr);
	
	pwtObj->mFilterBankArr=mFilterBankArr;
	pwtObj->freBandArr=freBandArr;
	pwtObj->binBandArr=binBandArr;
}

float *pwtObj_getFreBandArr(PWTObj pwtObj){

	return pwtObj->freBandArr;
}

int *pwtObj_getBinBandArr(PWTObj pwtObj){

	return pwtObj->binBandArr;
}

void pwtObj_pwt(PWTObj pwtObj,float *dataArr,float *mRealArr4,float *mImageArr4){

	__pwtObj_pwt(pwtObj,dataArr,pwtObj->mFilterBankArr,0,mRealArr4,mImageArr4);
}

// dataArr is NULL,use cwt middle result;
void pwtObj_pwtDet(PWTObj pwtObj,float *dataArr,float *mRealArr4,float *mImageArr4){

	if(pwtObj->mFilterBankDetArr){
		__pwtObj_pwt(pwtObj,dataArr,pwtObj->mFilterBankDetArr,1,mRealArr4,mImageArr4);
	}
}

void pwtObj_enableDet(PWTObj pwtObj,int flag){
	int fftLength=0; 
	int num=0;

	float *mFilterBankArr=NULL; // num*fftLength
	float *mFilterBankDetArr=NULL; 

	float *wArr=NULL;

	if(!flag||pwtObj->mFilterBankDetArr){
		return;
	}

	fftLength=pwtObj->fftLength;
	num=pwtObj->num;

	mFilterBankArr=pwtObj->mFilterBankArr;

	mFilterBankDetArr=__vnew(num*fftLength, NULL);
	wArr=__vnew(fftLength, NULL);
	for(int i=0;i<=fftLength/2;i++){
		wArr[i]=i*2*M_PI/fftLength;
	}
	for(int i=fftLength/2+1,j=fftLength/2-1;i<fftLength&&j>=0;i++,j--){
		wArr[i]=-wArr[j];
	}

	__mmul_vector(mFilterBankArr, wArr, 1, num, fftLength, 1, mFilterBankDetArr);

	pwtObj->mFilterBankDetArr=mFilterBankDetArr;

	// debug
	// {
	// 	printf("bankArr is :\n");
	// 	__mdebug(cwtObj->mFilterBankArr,  num, fftLength, 1);
	// 	printf("\n");

	// 	printf("det is :\n");
	// 	__mdebug(cwtObj->mFilterBankDetArr,  num, fftLength, 1);
	// 	printf("\n");
	// }

	free(wArr);
}

static void __pwtObj_pwt(PWTObj pwtObj,float *dataArr,float *mFilterBankArr,int iFlag,float *mRealArr4,float *mImageArr4){
	FFTObj fftObj=NULL;
	DFTObj dftObj=NULL;

	int fftLength=0; // 1<<radix2Exp+2*padLength
	int dataLength=0; // data length=1<<radix2Exp
	int padLength=0; // 0||dataLength/2||log2(dataLength) ???

	int num=0;

	float *realArr1=NULL; // fft data result
	float *imageArr1=NULL;

	float *mRealArr2=NULL; // mFilterBankArr.dot(fftData)  
	float *mImageArr2=NULL;

	float *mRealArr3=NULL; // ifft(mArr2) !isPad not use
	float *mImageArr3=NULL;

	float *curDataArr=NULL; 

	fftObj=pwtObj->fftObj;
	dftObj=pwtObj->dftObj;

	fftLength=pwtObj->fftLength;
	dataLength=pwtObj->dataLength;
	padLength=pwtObj->padLength;

	num=pwtObj->num;

	realArr1=pwtObj->realArr1;
	imageArr1=pwtObj->imageArr1;

	mRealArr2=pwtObj->mRealArr2;
	mImageArr2=pwtObj->mImageArr2;

	mRealArr3=pwtObj->mRealArr3;
	mImageArr3=pwtObj->mImageArr3;

	curDataArr=pwtObj->curDataArr;

	if(dataArr){
		// 1. deal padding data
		if(padLength){
			for(int i=padLength-1,j=0;i>=0;i--,j++){
				curDataArr[j]=dataArr[i];
			}

			memcpy(curDataArr+padLength, dataArr, sizeof(float )*dataLength);

			for(int i=dataLength-1,j=padLength+dataLength;i>=dataLength-padLength;i--,j++){
				curDataArr[j]=dataArr[i];
			}
		}

	 	// 2. fft
		if(fftObj){
			fftObj_fft(fftObj, padLength?curDataArr:dataArr, NULL, realArr1, imageArr1);
		}
		else{
			dftObj_dft(dftObj, padLength?curDataArr:dataArr, NULL, realArr1, imageArr1);
		}
	}
	
 	// 3. mFilterBankArr.dot(fftData) 
	for(int i=0;i<num;i++){
		for(int j=0;j<fftLength;j++){
			if(!iFlag){ // mFilterBankArr
				mRealArr2[j+i*fftLength]=mFilterBankArr[j+i*fftLength]*realArr1[j];
				mImageArr2[j+i*fftLength]=mFilterBankArr[j+i*fftLength]*imageArr1[j];
			}
			else{ // mFilterBankDetArr j
				mRealArr2[j+i*fftLength]=-mFilterBankArr[j+i*fftLength]*imageArr1[j];
				mImageArr2[j+i*fftLength]=mFilterBankArr[j+i*fftLength]*realArr1[j];
			}
		}
	}

 	// 4. ifft
	if(fftObj){
		if(padLength){ // has padding
			for(int i=0;i<num;i++){
				fftObj_ifft(fftObj, mRealArr2+i*fftLength, mImageArr2+i*fftLength, 
							mRealArr3+i*fftLength, mImageArr3+i*fftLength);
			}

			for(int i=0;i<num;i++){
				for(int j=padLength,k=0;j<padLength+dataLength;j++,k++){
					mRealArr4[k+i*dataLength]=mRealArr3[j+i*fftLength];
					mImageArr4[k+i*dataLength]=mImageArr3[j+i*fftLength];
				}
			}
		}
		else{
			for(int i=0;i<num;i++){
				fftObj_ifft(fftObj, mRealArr2+i*fftLength, mImageArr2+i*fftLength, 
							mRealArr4+i*fftLength, mImageArr4+i*fftLength);
			}
		}
	}
	else{
		if(padLength){ // has padding
			for(int i=0;i<num;i++){
				dftObj_idft(dftObj, mRealArr2+i*fftLength, mImageArr2+i*fftLength, 
							mRealArr3+i*fftLength, mImageArr3+i*fftLength);
			}

			for(int i=0;i<num;i++){
				for(int j=padLength,k=0;j<padLength+dataLength;j++,k++){
					mRealArr4[k+i*dataLength]=mRealArr3[j+i*fftLength];
					mImageArr4[k+i*dataLength]=mImageArr3[j+i*fftLength];
				}
			}
		}
		else{
			for(int i=0;i<num;i++){
				dftObj_idft(dftObj, mRealArr2+i*fftLength, mImageArr2+i*fftLength, 
							mRealArr4+i*fftLength, mImageArr4+i*fftLength);
			}
		}
	}
}

void pwtObj_free(PWTObj pwtObj){
	FFTObj fftObj=NULL;
	DFTObj dftObj=NULL;

	int *binBandArr=NULL;
	float *freBandArr=NULL;

	float *mFilterBankArr=NULL; // num*fftLength
	float *mFilterBankDetArr=NULL;

	float *realArr1=NULL; // fft data result
	float *imageArr1=NULL;

	float *mRealArr2=NULL; // mFilterBankArr.dot(fftData)  
	float *mImageArr2=NULL;

	float *mRealArr3=NULL; // ifft(mArr2) !isPad not use
	float *mImageArr3=NULL;

	float *curDataArr=NULL; 

	if(pwtObj){
		fftObj=pwtObj->fftObj;
		dftObj=pwtObj->dftObj;

		binBandArr=pwtObj->binBandArr;
		freBandArr=pwtObj->freBandArr;

		mFilterBankArr=pwtObj->mFilterBankArr;
		mFilterBankDetArr=pwtObj->mFilterBankDetArr;

		realArr1=pwtObj->realArr1;
		imageArr1=pwtObj->imageArr1;

		mRealArr2=pwtObj->mRealArr2;
		mImageArr2=pwtObj->mImageArr2;

		mRealArr3=pwtObj->mRealArr3;
		mImageArr3=pwtObj->mImageArr3;

		curDataArr=pwtObj->curDataArr;

		fftObj_free(fftObj);
		dftObj_free(dftObj);

		free(binBandArr);
		free(freBandArr);

		free(mFilterBankArr);
		free(mFilterBankDetArr);

		free(realArr1);
		free(imageArr1);

		free(mRealArr2);
		free(mImageArr2);

		free(mRealArr3);
		free(mImageArr3);

		free(curDataArr);

		free(pwtObj);
	}
}










