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
#include "filterbank/cwt_filterBank.h"

#include "cwt_algorithm.h"

struct OpaqueCWT{
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

	int binPerOctave; // for log

	float gamma;
	float beta;

	WaveletContinueType waveletType;
	SpectralFilterBankScaleType scaleType;

};

static void __cwtObj_init(CWTObj cwtObj);
static void __cwtObj_cwt(CWTObj cwtObj,float *dataArr,float *mFilterBankArr,int iFlag,float *mRealArr4,float *mImageArr4);

/***
	waveletType 
		'morse' gamma->wc+ beta->bin+ default
		'morlet' gamma->wc+ beta->bin-
		'bump' gamma->wc+ beta->bin-
		morlet->morse->bump ,fre resolution high,time resolution low
****/
int cwtObj_new(CWTObj *cwtObj,int num,int radix2Exp,
			 int *samplate,float *lowFre,float *highFre,int *binPerOctave,
			 WaveletContinueType *waveletType,
			 SpectralFilterBankScaleType *scaleType,
			 float *gamma,float *beta,
			 int *isPad){
	int status=0;
	CWTObj cwt=NULL;

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

	int _isPad=0;
	float _gamma=3;
	float _beta=20;

	int _flag=0;

	WaveletContinueType _waveletType=WaveletContinue_Morse;
	SpectralFilterBankScaleType _scaleType=SpectralFilterBankScale_Octave;

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

	if(waveletType){
		_waveletType=*waveletType;
	}
	
	if(scaleType){
		_scaleType=*scaleType;
		if(_scaleType>SpectralFilterBankScale_Log){
			printf("scaleType is error!\n");
			return 1;
		}
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

	if(_waveletType==WaveletContinue_Morlet){
		_gamma=6;
		_beta=2;
	}
	else if(_waveletType==WaveletContinue_Bump){
		_gamma=5;
		_beta=0.6;
	}
	else if(_waveletType==WaveletContinue_Paul){
		_gamma=4;
	}
	else if(_waveletType==WaveletContinue_DOG){
		_gamma=2;
		_beta=2;
	}
	else if(_waveletType==WaveletContinue_Mexican){
		_beta=2;
	}
	else if(_waveletType==WaveletContinue_Hermit){
		_gamma=5;
		_beta=2;
	}
	else if(_waveletType==WaveletContinue_Ricker){
		_gamma=4;
	}

	if(gamma){
		if(*gamma>0){
			_gamma=*gamma;
			
			if(_waveletType==WaveletContinue_DOG){
				int p=0;

				p=roundf(_gamma);
				if(p%2==0){
					_gamma=p;
				}
				else{
					_gamma=2;
				}
			}
		}
	}

	if(beta){
		if(*beta>0){
			_beta=*beta;
		}
	}

	if(isPad){
		_isPad=*isPad;
	}

	if(num<2||num>fftLength/2+1){ 
		printf("num is error!\n");
		return -1;
	}

	cwt=*cwtObj=(CWTObj )calloc(1, sizeof(struct OpaqueCWT ));

	dataLength=1<<radix2Exp;
	if(_isPad){
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

	cwt->fftObj=fftObj;
	cwt->dftObj=dftObj;

	cwt->fftLength=fftLength;
	cwt->dataLength=dataLength;
	cwt->padLength=padLength;

	cwt->num=num;

	cwt->realArr1=realArr1;
	cwt->imageArr1=imageArr1;

	cwt->mRealArr2=mRealArr2;
	cwt->mImageArr2=mImageArr2;

	cwt->mRealArr3=mRealArr3;
	cwt->mImageArr3=mImageArr3;

	cwt->fftObj=fftObj;

	cwt->curDataArr=curDataArr;

	cwt->samplate=_samplate;
	cwt->lowFre=_lowFre;
	cwt->highFre=_highFre;

	cwt->binPerOctave=_binPerOctave;

	cwt->gamma=_gamma;
	cwt->beta=_beta;

	cwt->waveletType=_waveletType;
	cwt->scaleType=_scaleType;

	__cwtObj_init(cwt);

	return status;
}

float *cwtObj_getFreBandArr(CWTObj cwtObj){

	return cwtObj->freBandArr;
}

int *cwtObj_getBinBandArr(CWTObj cwtObj){

	return cwtObj->binBandArr;
}

void cwtObj_cwt(CWTObj cwtObj,float *dataArr,float *mRealArr4,float *mImageArr4){
	
	__cwtObj_cwt(cwtObj,dataArr,cwtObj->mFilterBankArr,0,mRealArr4,mImageArr4);

}

// dataArr is NULL,use cwt middle result;
void cwtObj_cwtDet(CWTObj cwtObj,float *dataArr,float *mRealArr4,float *mImageArr4){

	if(cwtObj->mFilterBankDetArr){
		__cwtObj_cwt(cwtObj,dataArr,cwtObj->mFilterBankDetArr,1,mRealArr4,mImageArr4);
	}
}

// mFilterBankArr&mFilterBankDetArr
static void __cwtObj_cwt(CWTObj cwtObj,float *dataArr,float *mFilterBankArr,int iFlag,float *mRealArr4,float *mImageArr4){
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

	fftObj=cwtObj->fftObj;
	dftObj=cwtObj->dftObj;

	fftLength=cwtObj->fftLength;
	dataLength=cwtObj->dataLength;
	padLength=cwtObj->padLength;

	num=cwtObj->num;

	realArr1=cwtObj->realArr1;
	imageArr1=cwtObj->imageArr1;

	mRealArr2=cwtObj->mRealArr2;
	mImageArr2=cwtObj->mImageArr2;

	mRealArr3=cwtObj->mRealArr3;
	mImageArr3=cwtObj->mImageArr3;

	curDataArr=cwtObj->curDataArr;

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

void cwtObj_enableDet(CWTObj cwtObj,int flag){
	int fftLength=0; 
	int num=0;

	float *mFilterBankArr=NULL; // num*fftLength
	float *mFilterBankDetArr=NULL; 

	float *wArr=NULL;

	if(!flag||cwtObj->mFilterBankDetArr){
		return;
	}

	fftLength=cwtObj->fftLength;
	num=cwtObj->num;

	mFilterBankArr=cwtObj->mFilterBankArr;

	mFilterBankDetArr=__vnew(num*fftLength, NULL);
	wArr=__vnew(fftLength, NULL);
	for(int i=0;i<=fftLength/2;i++){
		wArr[i]=i*2*M_PI/fftLength;
	}
	for(int i=fftLength/2+1,j=fftLength/2-1;i<fftLength&&j>=0;i++,j--){
		wArr[i]=-wArr[j];
	}

	__mmul_vector(mFilterBankArr, wArr, 1, num, fftLength, 1, mFilterBankDetArr);

	cwtObj->mFilterBankDetArr=mFilterBankDetArr;

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

void cwtObj_free(CWTObj cwtObj){
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

	if(cwtObj){
		fftObj=cwtObj->fftObj;
		dftObj=cwtObj->dftObj;

		binBandArr=cwtObj->binBandArr;
		freBandArr=cwtObj->freBandArr;

		mFilterBankArr=cwtObj->mFilterBankArr;
		mFilterBankDetArr=cwtObj->mFilterBankDetArr;

		realArr1=cwtObj->realArr1;
		imageArr1=cwtObj->imageArr1;

		mRealArr2=cwtObj->mRealArr2;
		mImageArr2=cwtObj->mImageArr2;

		mRealArr3=cwtObj->mRealArr3;
		mImageArr3=cwtObj->mImageArr3;

		curDataArr=cwtObj->curDataArr;

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

		free(cwtObj);
	}
}

static void __cwtObj_init(CWTObj cwtObj){
	int *binBandArr=NULL;
	float *freBandArr=NULL;

	float *mFilterBankArr=NULL;

	int num=0;
	int dataLength=0;
	int samplate=0;
	int padLength=0;
	int fftLength=0;

	WaveletContinueType waveletType;
	float gamma=0;
	float beta=0;

	SpectralFilterBankScaleType scaleType;

	float lowFre=0;
	float highFre=0;
	int binPerOctave=0;

	num=cwtObj->num;
	dataLength=cwtObj->dataLength;
	samplate=cwtObj->samplate;
	padLength=cwtObj->padLength;
	fftLength=cwtObj->fftLength;

	waveletType=cwtObj->waveletType;
	gamma=cwtObj->gamma;
	beta=cwtObj->beta;

	scaleType=cwtObj->scaleType;

	lowFre=cwtObj->lowFre;
	highFre=cwtObj->highFre;
	binPerOctave=cwtObj->binPerOctave;

	mFilterBankArr=__vnew(num*fftLength, NULL);
	freBandArr=__vnew(num+2, NULL);
	binBandArr=__vnewi(num+2, NULL);

	cwt_filterBank(num,dataLength,samplate,padLength,
				waveletType,gamma,beta,
				scaleType,
				lowFre,highFre,binPerOctave,
				mFilterBankArr,
				freBandArr,
				binBandArr);

	cwtObj->mFilterBankArr=mFilterBankArr;
	cwtObj->freBandArr=freBandArr;
	cwtObj->binBandArr=binBandArr;

	// debug
	// {
	// 	printf("filterBank is :\n");
	// 	__mdebug(mFilterBankArr, num, fftLength, 1);
	// 	printf("\n");
	// }
}









