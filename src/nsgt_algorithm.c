// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "dsp/flux_window.h"
#include "dsp/fft_algorithm.h"
#include "dsp/dft_algorithm.h"

#include "filterbank/auditory_filterBank.h"
#include "filterbank/nsgt_filterBank.h"

#include "nsgt_algorithm.h"

struct OpaqueNSGT{
	FFTObj fftObj;
	int fftLength; // data length

	int num;
	float *windowDataArr;
	int *windowLengthArr; // timeLengthArr

	int *binBandArr;
	float *freBandArr;
	int *offsetArr;

	DFTObj *dftArr;
	int *dftLenArr;
	int dftLength;

	int maxWindowLength;
	int totalWindowLength;
	int minWindowLength;

	float *realArr1; // fft data result
	float *imageArr1;

	float *realArr2; // reassign freData*windowData =>n/2~n,1~n/2
	float *imageArr2;

	float *realArr3; // cell result data
	float *imageArr3;

	float *maxTimeArr; // maxWindowLength+1 
	float **timeArrArr;

	// params
	int samplate; // filterBank 
	float lowFre;
	float highFre;

	int binPerOctave; // for Log/LogChroma

	NSGTFilterBankType nsgtFilterBankType;
	SpectralFilterBankScaleType filterScaleType;
	SpectralFilterBankStyleType filterStyleType;
	SpectralFilterBankNormalType filterNormalType;

};

static void __nsgtObj_init(NSGTObj nsgtObj);
static void __nsgtObj_dealTime(NSGTObj nsgtObj);
static void __nsgtObj_dealDFT(NSGTObj nsgtObj);

static int __arr_has(int *arr,int length,int value);
static int __arr_getIndex(int *arr,int length,int value);

int nsgtObj_new(NSGTObj *nsgtObj,int num,int radix2Exp,
				int *samplate,float *lowFre,float *highFre,int *binPerOctave,
				int *minLength,
				NSGTFilterBankType *nsgtFilterBankType,
				SpectralFilterBankScaleType *filterScaleType,
				SpectralFilterBankStyleType *filterStyleType,
				SpectralFilterBankNormalType *filterNormalType){
	int status=0;
	NSGTObj nsgt=NULL;

	FFTObj fftObj=NULL;
	int fftLength=0;

	float *realArr1=NULL; // fft data result
	float *imageArr1=NULL;

	float *realArr2=NULL; // reassign freData*windowData =>n/2~n,1~n/2
	float *imageArr2=NULL;

	int _samplate=32000; // filterBank 
	float _lowFre=0;
	float _highFre=0;

	int _binPerOctave=12; // log

	int minWindowLength=3;

	NSGTFilterBankType _nsgtFilterBankType=NSGTFilterBank_Efficient;
	SpectralFilterBankScaleType _filterScaleType=SpectralFilterBankScale_Octave;
	SpectralFilterBankStyleType _filterStyleType=SpectralFilterBankStyle_Hann;
	SpectralFilterBankNormalType _filterNormalType=SpectralFilterBankNormal_BandWidth;

	if(minLength){
		if(*minLength>0){
			minWindowLength=*minLength;
		}
	}

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

	if(nsgtFilterBankType){
		_nsgtFilterBankType=*nsgtFilterBankType;
	}

	if(filterScaleType){
		_filterScaleType=*filterScaleType;
		if(_filterScaleType>SpectralFilterBankScale_Log){
			printf("scaleType is error!\n");
			return 1;
		}
	}

	if(filterStyleType){
		_filterStyleType=*filterStyleType;
		if(_filterStyleType==SpectralFilterBankStyle_Gammatone){
			_filterStyleType=SpectralFilterBankStyle_Hann;
		}
	}
	
	if(filterNormalType){
		_filterNormalType=*filterNormalType;
		if(_filterNormalType==SpectralFilterBankNormal_Area){
			 _filterNormalType=SpectralFilterBankNormal_BandWidth;
		}
	}

	_highFre=_samplate/2.0;

	if(lowFre){
		if(*lowFre>=0&&*lowFre<_samplate/2.0){
			_lowFre=*lowFre;
		}
	}

	if(_lowFre==0){
		if(_filterScaleType==SpectralFilterBankScale_Octave||
			_filterScaleType==SpectralFilterBankScale_Log){ // Log/Logspace

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
		if(_filterScaleType==SpectralFilterBankScale_Octave||
			_filterScaleType==SpectralFilterBankScale_Log){ // Log/Logspace
			
			_lowFre=powf(2, -45/12.0)*440;
			_highFre=powf(2, 38/12.0)*440;
		}
	}

	if(binPerOctave){
		if(*binPerOctave>=4&&*binPerOctave<=48){
			_binPerOctave=*binPerOctave;
		}
	}

	if(_filterScaleType==SpectralFilterBankScale_Linear){ // linear
		float det=0;

		det=_samplate/(float )fftLength;
		auditory_reviseLinearFre(num, _lowFre,_highFre, det,1,&_lowFre, &_highFre);

		if(_highFre>_samplate/2.0){
			printf("scale linear: lowFre and num is large, overflow error\n");
			return -1;
		}
	}
	else if(_filterScaleType==SpectralFilterBankScale_Octave){ // log
		auditory_reviseLogFre(num, _lowFre,_highFre, _binPerOctave,1, &_lowFre, &_highFre);

		if(_highFre>_samplate/2.0){
			printf("scale log: lowFre and num is large, overflow error!\n");
			return -1;
		}
	}

	if(num<2||num>fftLength/2+1){ 
		printf("num is error!\n");
		return -1;
	}

	realArr1=__vnew(fftLength, NULL);
	imageArr1=__vnew(fftLength, NULL);

	realArr2=__vnew(fftLength, NULL);
	imageArr2=__vnew(fftLength, NULL);

	nsgt=*nsgtObj=(NSGTObj )calloc(1, sizeof(struct OpaqueNSGT ));
	fftObj_new(&fftObj, radix2Exp);

	nsgt->minWindowLength=minWindowLength;

	nsgt->fftObj=fftObj;
	nsgt->fftLength=fftLength;
	nsgt->num=num;

	nsgt->realArr1=realArr1;
	nsgt->imageArr1=imageArr1;

	nsgt->realArr2=realArr2;
	nsgt->imageArr2=imageArr2;

	nsgt->samplate=_samplate;
	nsgt->lowFre=_lowFre;
	nsgt->highFre=_highFre;
	nsgt->binPerOctave=_binPerOctave;

	nsgt->nsgtFilterBankType=_nsgtFilterBankType;
	nsgt->filterScaleType=_filterScaleType;
	nsgt->filterStyleType=_filterStyleType;
	nsgt->filterNormalType=_filterNormalType;

	__nsgtObj_init(nsgt);
	__nsgtObj_dealTime(nsgt);

	return status;
}

static void __nsgtObj_dealTime(NSGTObj nsgtObj){
	int num=0;
	int *windowLengthArr=NULL;

	int maxWindowLength=0;
	int fftLength=0; 
	int samplate=0;

	float *maxTimeArr=NULL; // maxWindowLength+1 
	float **timeArrArr=NULL;

	float time=0;
	float curLen=0;

	num=nsgtObj->num;
	windowLengthArr=nsgtObj->windowLengthArr;

	maxWindowLength=nsgtObj->maxWindowLength;
	fftLength=nsgtObj->fftLength;
	samplate=nsgtObj->samplate;

	time=fftLength/(float )samplate;
	maxTimeArr=__vlinspace(0, time, maxWindowLength+1, 0);
	timeArrArr=(float **)calloc(num, sizeof(float *));
	for(int i=0;i<num;i++){
		float det=0;
		float offset=0;

		curLen=windowLengthArr[i];
		det=(curLen-2>=0?curLen-2:0);
		offset=time/(curLen+det);

		timeArrArr[i]=__vlinspace(-offset, time+offset, curLen+1, 0);
	}

	nsgtObj->maxTimeArr=maxTimeArr;
	nsgtObj->timeArrArr=timeArrArr;
}

static void __nsgtObj_dealDFT(NSGTObj nsgtObj){
	DFTObj *dftArr=NULL;
	int *dftLenArr=NULL;
	int dftLength=0;

	int num=0;
	int *windowLengthArr=NULL; 

	int flag=0;

	num=nsgtObj->num;
	windowLengthArr=nsgtObj->windowLengthArr;

	dftArr=(DFTObj *)calloc(num, sizeof(DFTObj ));
	dftLenArr=__vnewi(num, NULL);
	for(int i=0;i<num;i++){
		flag=__arr_has(dftLenArr, dftLength, windowLengthArr[i]);
		if(!flag){
			dftLenArr[dftLength]=windowLengthArr[i];
			dftLength++;
		}
	}

	for(int i=0;i<dftLength;i++){
		dftObj_new(dftArr+i, dftLenArr[i]);
	}

	nsgtObj->dftArr=dftArr;
	nsgtObj->dftLenArr=dftLenArr;
	nsgtObj->dftLength=dftLength;
}

static void __nsgtObj_init(NSGTObj nsgtObj){
	int fftLength=0; // data length
	int num=0;

	int samplate=0; // filterBank 
	float lowFre=0;
	float highFre=0;
	int binPerOctave=0; // for Log/LogChroma

	int minWindowLength=0;

	NSGTFilterBankType nsgtFilterBankType;
	SpectralFilterBankScaleType filterScaleType;
	SpectralFilterBankStyleType filterStyleType;
	SpectralFilterBankNormalType filterNormalType;

	float *windowDataArr=NULL;
	int *windowLengthArr=NULL; // timeLengthArr

	int *binBandArr=NULL;
	float *freBandArr=NULL;
	int *offsetArr=NULL;

	int totalWindowLength=0;
	int maxWindowLength=0; 

	float *realArr3=NULL; // cell data result
	float *imageArr3=NULL;

	int type=0;

	minWindowLength=nsgtObj->minWindowLength;

	fftLength=nsgtObj->fftLength;
	num=nsgtObj->num;

	samplate=nsgtObj->samplate;
	lowFre=nsgtObj->lowFre;
	highFre=nsgtObj->highFre;
	binPerOctave=nsgtObj->binPerOctave;

	nsgtFilterBankType=nsgtObj->nsgtFilterBankType;
	filterScaleType=nsgtObj->filterScaleType;
	filterStyleType=nsgtObj->filterStyleType;
	filterNormalType=nsgtObj->filterNormalType;

	if(nsgtFilterBankType==NSGTFilterBank_Standard){
		type=1;
	}

	windowLengthArr=__vnewi(num, NULL);
	binBandArr=__vnewi(num, NULL);
	freBandArr=__vnew(num, NULL);
	offsetArr=__vnewi(num, NULL);

	nsgt_filterBank(num,fftLength,samplate, minWindowLength,
					type,
					filterScaleType,
					filterStyleType,
					filterNormalType,
					lowFre,highFre,binPerOctave,
					&windowDataArr,windowLengthArr,
					freBandArr,binBandArr,offsetArr,
					&maxWindowLength,&totalWindowLength);

	realArr3=__vnew(totalWindowLength, NULL);
	imageArr3=__vnew(totalWindowLength, NULL);

	nsgtObj->windowDataArr=windowDataArr;
	nsgtObj->windowLengthArr=windowLengthArr;

	nsgtObj->binBandArr=binBandArr;
	nsgtObj->freBandArr=freBandArr;
	nsgtObj->offsetArr=offsetArr;

	nsgtObj->maxWindowLength=maxWindowLength;
	nsgtObj->totalWindowLength=totalWindowLength;

	nsgtObj->realArr3=realArr3;
	nsgtObj->imageArr3=imageArr3;

	__nsgtObj_dealDFT(nsgtObj);

	// debug
	// {
	// 	printf("windowDataArr is \n");
	// 	int index=0;
	// 	for(int i=0;i<num;i++){
	// 		int len=0;

	// 		len=windowLengthArr[i];
	// 		printf(" %d: %d\n",i,len);
	// 		printf("	");
	// 		for(int j=0;j<len;j++){
	// 			printf("%d: %f, ",j,windowDataArr[index+j]);
	// 		}
	// 		printf("\n\n");

	// 		index+=len;
	// 	}
	// }

}

// default 3  minLength>=1
void nsgtObj_setMinLength(NSGTObj nsgtObj,int minLength){
	int minWindowLength=0;

	float *windowDataArr=NULL;
	int *windowLengthArr=NULL; // timeLengthArr

	int *binBandArr=NULL;
	float *freBandArr=NULL;
	int *offsetArr=NULL;

	float *realArr3=NULL; // cell data result
	float *imageArr3=NULL;

	DFTObj *dftArr=NULL;
	int *dftLenArr=NULL;
	int dftLength=0;

	minWindowLength=nsgtObj->minWindowLength;
	if(minLength!=minWindowLength&&minLength>=1){
		windowDataArr=nsgtObj->windowDataArr;
		windowLengthArr=nsgtObj->windowLengthArr;

		binBandArr=nsgtObj->binBandArr;
		freBandArr=nsgtObj->freBandArr;
		offsetArr=nsgtObj->offsetArr;

		realArr3=nsgtObj->realArr3;
		imageArr3=nsgtObj->imageArr3;

		dftArr=nsgtObj->dftArr;
		dftLenArr=nsgtObj->dftLenArr;
		dftLength=nsgtObj->dftLength;

		free(windowDataArr);
		free(windowLengthArr);

		free(binBandArr);
		free(freBandArr);
		free(offsetArr);

		free(realArr3);
		free(imageArr3);

		for(int i=0;i<dftLength;i++){
			dftObj_free(dftArr[i]);
		}
		free(dftArr);
		free(dftLenArr);

		nsgtObj->minWindowLength=minLength;
		__nsgtObj_init(nsgtObj);
	}
}

void nsgtObj_nsgt(NSGTObj nsgtObj,float *dataArr,float *mRealArr3,float *mImageArr3){
	FFTObj fftObj=NULL;
	int fftLength=0;

	int num=0;
	float *windowDataArr=NULL;
	int *windowLengthArr=NULL; 

	float *realArr1=NULL; // data fft/dft result
	float *imageArr1=NULL;

	float *realArr2=NULL; // data fft/dft result
	float *imageArr2=NULL;

	float *realArr3=NULL; // data fft/dft result
	float *imageArr3=NULL;

	int *binBandArr=NULL;
	int *offsetArr=NULL;

	int maxWindowLength=0;
	float *maxTimeArr=NULL;
	float **timeArrArr=NULL;

	DFTObj *dftArr=NULL;
	int *dftLenArr=NULL;
	int dftLength=0;

	int curLen=0;
	int index=0;

	fftObj=nsgtObj->fftObj;
	fftLength=nsgtObj->fftLength;

	num=nsgtObj->num;
	windowDataArr=nsgtObj->windowDataArr;
	windowLengthArr=nsgtObj->windowLengthArr;

	realArr1=nsgtObj->realArr1;
	imageArr1=nsgtObj->imageArr1;

	realArr2=nsgtObj->realArr2;
	imageArr2=nsgtObj->imageArr2;

	realArr3=nsgtObj->realArr3;
	imageArr3=nsgtObj->imageArr3;

	binBandArr=nsgtObj->binBandArr;
	offsetArr=nsgtObj->offsetArr;

	maxWindowLength=nsgtObj->maxWindowLength;
	maxTimeArr=nsgtObj->maxTimeArr;
	timeArrArr=nsgtObj->timeArrArr;
	
	dftArr=nsgtObj->dftArr;
	dftLenArr=nsgtObj->dftLenArr;
	dftLength=nsgtObj->dftLength;

	// 1. fft(dataLength)
	fftObj_fft(fftObj, dataArr, NULL, realArr1, imageArr1);

	// 2. ifft
	for(int i=0;i<num;i++){
		int _dIndex=0;
		int _offset=0;
		float _value=0;


		curLen=windowLengthArr[i];
		_offset=offsetArr[i];
		for(int j=0,k=curLen-curLen/2;j<curLen;j++,k++){
			// float r1=0;
			// float i1=0;

			// if(_offset>=0&&_offset<fftLength){
			// 	r1=realArr1[_offset];
			// 	i1=imageArr1[_offset];
			// }

			_value=windowDataArr[index+j];
			if(k>=curLen){
				k=0;
			}

			if(_offset>fftLength-1){
				_offset=fftLength-1;
			}
			else if(_offset<0){
				_offset=0;
			}

			realArr2[k]=realArr1[_offset]*_value;
			imageArr2[k]=imageArr1[_offset]*_value;

			_offset++;
		}

		_dIndex=__arr_getIndex(dftLenArr, dftLength, curLen);
		dftObj_idft(dftArr[_dIndex], realArr2, imageArr2, realArr3+index, imageArr3+index);
		index+=curLen;
	}

	// 3. to matrix
	index=0;
	for(int i=0;i<num;i++){
		int start=0;

		curLen=windowLengthArr[i];
		for(int j=0;j<maxWindowLength;j++){
			for(int k=start;k<curLen+1;k++){
				if(maxTimeArr[j]<timeArrArr[i][k]){
					mRealArr3[i*maxWindowLength+j]=realArr3[index+k-1];
					mImageArr3[i*maxWindowLength+j]=imageArr3[index+k-1];
					
					start=k;
					break;
				}
			}
		}
				
		index+=curLen;
	}
}

// test cell data
void nsgtObj_getCellData(NSGTObj nsgtObj,float **realArr3,float **imageArr3){

	*realArr3=nsgtObj->realArr3;
	*imageArr3=nsgtObj->imageArr3;
}

int nsgtObj_getMaxTimeLength(NSGTObj nsgtObj){

	return nsgtObj->maxWindowLength;
}

int nsgtObj_getTotalTimeLength(NSGTObj nsgtObj){

	return nsgtObj->totalWindowLength;
}

int *nsgtObj_getTimeLengthArr(NSGTObj nsgtObj){

	return nsgtObj->windowLengthArr;
}

float *nsgtObj_getFreBandArr(NSGTObj nsgtObj){

	return nsgtObj->freBandArr;
}

int *nsgtObj_getBinBandArr(NSGTObj nsgtObj){

	return nsgtObj->binBandArr;
}

void nsgtObj_free(NSGTObj nsgtObj){
	FFTObj fftObj=NULL;

	int num=0;
	float *windowDataArr=NULL;
	int *windowLengthArr=NULL; // timeLengthArr

	int *binBandArr=NULL;
	float *freBandArr=NULL;
	int *offsetArr=NULL;

	DFTObj *dftArr=NULL;
	int *dftLenArr=NULL;
	int dftLength=0;

	float *realArr1=NULL; // fft data result
	float *imageArr1=NULL;

	float *realArr2=NULL; // reassign freData*windowData =>n/2~n,1~n/2
	float *imageArr2=NULL;

	float *realArr3=NULL; // cell result data
	float *imageArr3=NULL;

	float *maxTimeArr=NULL; // maxWindowLength+1 
	float **timeArrArr=NULL;

	if(!nsgtObj){
		return;
	}

	fftObj=nsgtObj->fftObj;

	num=nsgtObj->num;
	windowDataArr=nsgtObj->windowDataArr;
	windowLengthArr=nsgtObj->windowLengthArr;

	binBandArr=nsgtObj->binBandArr;
	freBandArr=nsgtObj->freBandArr;
	offsetArr=nsgtObj->offsetArr;

	dftArr=nsgtObj->dftArr;
	dftLenArr=nsgtObj->dftLenArr;
	dftLength=nsgtObj->dftLength;

	realArr1=nsgtObj->realArr1;
	imageArr1=nsgtObj->imageArr1;

	realArr2=nsgtObj->realArr2;
	imageArr2=nsgtObj->imageArr2;

	realArr3=nsgtObj->realArr3;
	imageArr3=nsgtObj->imageArr3;

	maxTimeArr=nsgtObj->maxTimeArr;
	timeArrArr=nsgtObj->timeArrArr;

	fftObj_free(fftObj);

	free(windowDataArr);
	free(windowLengthArr);

	free(binBandArr);
	free(freBandArr);
	free(offsetArr);

	for(int i=0;i<dftLength;i++){
		dftObj_free(dftArr[i]);
	}
	free(dftArr);
	free(dftLenArr);

	free(realArr1);
	free(imageArr1);

	free(realArr2);
	free(imageArr2);

	free(realArr3);
	free(imageArr3);

	free(maxTimeArr);
	for(int i=0;i<num;i++){
		free(timeArrArr[i]);
	}
	free(timeArrArr);

	free(nsgtObj);
}

static int __arr_has(int *arr,int length,int value){
	int flag=0;

	if(arr&&length){
		for(int i=0;i<length;i++){
			if(arr[i]==value){
				flag=1;
				break;
			}
		}
	}

	return flag;
}

static int __arr_getIndex(int *arr,int length,int value){
	int index=-1;

	for(int i=0;i<length;i++){
		if(arr[i]==value){
			index=i;
			break;
		}
	}

	return index;
}








