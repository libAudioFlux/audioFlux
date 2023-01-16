// clang 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"

#include "cqt_filterBank.h"

// tempKernel
static void __cqt_calTempArr(int num,float *freBandArr, int samplate,
				float *lenArr,int fftLength,
				SpectralFilterBankNormalType normType,WindowType *winType,
				float *mReaArr,float *mImageArr);

// 0 fftLength/2+1 1 fftLength
static void __cqt_filterBank(int num,float *freBandArr,int samplate,
				int binPerOctave,SpectralFilterBankNormalType normType,WindowType *winType,
				float *factor,float *beta,float *thresh,
				float *lenArr,int fftLength,
				float *mRealFilterBankArr,float *mImageFilterBankArr,int type);

// filterBank相关
/***
	cqt filterBank 基于音乐octave分割的freBandArr
		normType 归一化方式 ???
		winType 'hann' paper 'hamm'
		factor 1 调节时间分辨率 
		beta 0 调节可变Q
		thresh 0.01 paper 0.0054
	1. tempKernel计算
	2. freKernel计算
	3. 阈值过滤
	mFilterBankArr=>num*(fftLength/2+1) num即length 非稀疏矩阵 
****/
void cqt_filterBank(int num,float *freBandArr,int samplate,
				int binPerOctave,SpectralFilterBankNormalType normType,WindowType *winType,
				float *factor,float *beta,float *thresh,
				float *lenArr,int fftLength,
				float *mRealFilterBankArr,float *mImageFilterBankArr){

	__cqt_filterBank(num,freBandArr,samplate,
					binPerOctave,normType,winType,
					factor,beta,thresh,
					lenArr,fftLength,
					mRealFilterBankArr,mImageFilterBankArr,
					0);
}

void cqt_downFilterBank(int num,float *freBandArr,int samplate,
				int binPerOctave,SpectralFilterBankNormalType normType,WindowType *winType,
				float *factor,float *beta,float *thresh,
				float *lenArr,int fftLength,
				float *mRealFilterBankArr,float *mImageFilterBankArr){
	int octaveNum=0; // downNum

	FFTObj fftObj=NULL;
	int m=0;

	float *mRealTempArr=NULL; // tempKernel
	float *mImageTempArr=NULL;

	float *_mRealArr1=NULL;
	float *_mImageArr1=NULL;

	float _factor=1;
	float _beta=0;
	float _thresh=0.01; // paper 0.0054
	
	if(!mRealFilterBankArr||!mImageFilterBankArr){
		return;
	}

	if(factor){
		_factor=*factor;
	}

	if(beta){
		_beta=*beta;
	}

	if(thresh){
		_thresh=*thresh;
	}

	octaveNum=num/binPerOctave;

	_mRealArr1=__vnew(num*fftLength, NULL);
	_mImageArr1=__vnew(num*fftLength, NULL);

	m=util_powerTwoBit(fftLength);
	fftObj_new(&fftObj, m);

	mRealTempArr=__vnew(binPerOctave*fftLength, NULL);
	mImageTempArr=__vnew(binPerOctave*fftLength, NULL);
	for(int i=octaveNum-1;i>=0;i--){
		// 1. tempKernel计算
		__cqt_calTempArr(binPerOctave,freBandArr+i*binPerOctave,samplate,
					lenArr,fftLength,
					normType,winType,
					mRealTempArr,mImageTempArr);

		// 2. freKernel计算
		for(int j=0;j<binPerOctave;j++){
			int index1=0;
			int index2=0;

			index1=j*fftLength;
			index2=(i*binPerOctave+j)*fftLength;
			fftObj_fft(fftObj,
					mRealTempArr+index1,mImageTempArr+index1,
					_mRealArr1+index2,_mImageArr1+index2);
		}

		samplate/=2;
	}

	// 3.阈值过滤
	_thresh*=_thresh;
	for(int i=0;i<num;i++){
		for(int j=0;j<fftLength/2+1;j++){
			float _value1=0;
			float _value2=0;

			_value1=_mRealArr1[i*fftLength+j];
			_value2=_mImageArr1[i*fftLength+j];
			if(_value1*_value1+_value2*_value2>_thresh){
				mRealFilterBankArr[i*(fftLength/2+1)+j]=_value1;
				mImageFilterBankArr[i*(fftLength/2+1)+j]=_value2;
			}
		}
	}

	free(mRealTempArr);
	free(mImageTempArr);

	free(_mRealArr1);
	free(_mImageArr1);

	fftObj_free(fftObj);
}

// factor 1.0  sacle/(2^(1/binPerOctave)-1)
float cqt_calQ(int binPerOctave,float factor){
	float q=0;

	q=factor/(powf(2, 1.0/binPerOctave)-1);
	return q;
}

// num=octaveNum*binPerOctave;
float *cqt_calFreArr(float minFre,int num,int binPerOctave){
	float *arr=NULL;

	int octaveNum=0;

	float _fre=0;
	float _value=0;

	arr=__vnew(num+2, NULL);

	octaveNum=num/binPerOctave;
	_value=powf(2, 1.0/binPerOctave);

	// arr[0]=minFre/_value;
	for(int i=0;i<octaveNum;i++){
		_fre=minFre*(1<<i);
		arr[i*binPerOctave]=_fre;
		for(int j=1;j<binPerOctave;j++){
			_fre*=_value;
			arr[i*binPerOctave+j]=_fre;
		}
	}	
	// arr[num]=arr[num-1]*_value;

	return arr;
}

// ceil(Q*fs/freArr[i])
void cqt_calLengthArr(int num,float *freBandArr,int samplate,
				int binPerOctave,
				float *factor,float *beta,
				float *lenArr){
	float q=0;
	float value=0;

	float _factor=1;
	float _beta=0;

	if(factor){
		if(*factor>0){
			_factor=*factor;
		}
	}

	if(beta){
		_beta=*beta;
	}

	value=powf(2, 1.0/binPerOctave)-1;
	q=_factor/value;
	for(int i=0;i<num;i++){
		// lenArr[i]=ceilf(q*samplate/(freBandArr[i]+_beta/value));
		lenArr[i]=q*samplate/(freBandArr[i]+_beta/value);
	}
}

int cqt_calFFTLength(float minFre,int samplate,
				int binPerOctave,
				float *factor,float *beta){
	int fftLength=0;

	float q=0;
	float value=0;

	float _factor=1;
	float _beta=0;
	int _len=0;

	if(factor){
		if(*factor>0){
			_factor=*factor;
		}
	}

	if(beta){
		_beta=*beta;
	}

	value=powf(2, 1.0/binPerOctave)-1;
	q=_factor/value;

	_len=ceilf(q*samplate/(minFre+_beta/value));
	fftLength=util_ceilPowerTwo(_len);

	return fftLength;
}

/***
	1/Nk*W[Nk]*e^(2*PI*j*n*Q/Nk) n<Nk
	W[Nk] => winType 'hann'
	1/Nk => normType
	mRealArr => num*fftLength
	1. n-range 0~n;-n/2~n/2 2. norm 1/len;util.norm 3. padding zero/center
****/
static void __cqt_calTempArr(int num,float *freBandArr, int samplate,
				float *lenArr,int fftLength,
				SpectralFilterBankNormalType normType,WindowType *winType,
				float *mReaArr,float *mImageArr){
	WindowType _winType=Window_Hann;

	if(winType){
		if(*winType!=Window_Rect){
			_winType=*winType;
		}
	}

	for(int i=0;i<num;i++){
		int len=0;
		float fre=0;

		float *arr=NULL;
		float *wArr=NULL;

		int startIndex=0;
		float weight=0;
		
		len=ceilf(lenArr[i]);
		fre=freBandArr[i];

		__varange(0, len, 1, &arr); 
		// __varange(-len/2, len/2, 1, &arr);
		// wArr=window_createHann(len, 1); // ???
		wArr=window_calFFTWindow(_winType,len);

		weight=1;
		/***
			startIndex对应arr范围 0->0~len; ((fftLength-len)/2)->-en/2~len/2 ???
			应按公式标准0~len;然后居中两边填充
		****/
		startIndex=(fftLength-len)/2;
		for(int j=0;j<len;j++){
			float value=0;

			if(normType==SpectralFilterBankNormal_None){ // same hight
				weight=lenArr[i];
			}

			value=2*M_PI*arr[j]*fre/samplate;
			mReaArr[i*fftLength+j+startIndex]=cosf(value)*wArr[j]/weight; 
			mImageArr[i*fftLength+j+startIndex]=sinf(value)*wArr[j]/weight;
		}

		// norm => 统一类auditory 归一化处理
		weight=0; // reset
		if(normType==SpectralFilterBankNormal_Area){ // sum
			for(int j=0;j<len;j++){
				float value1=0;
				float value2=0;
				
				value1=mReaArr[i*fftLength+j+startIndex];
				value2=mImageArr[i*fftLength+j+startIndex];

				weight+=sqrtf(value1*value1+value2*value2);
			}

			for(int j=0;j<len;j++){
				mReaArr[i*fftLength+j+startIndex]/=weight; 
				mImageArr[i*fftLength+j+startIndex]/=weight;
			}
		}
		else if(normType==SpectralFilterBankNormal_BandWidth){
			weight=(freBandArr[i+1]-freBandArr[i-1])/2;
			for(int j=0;j<len;j++){
				mReaArr[i*fftLength+j+startIndex]/=weight; 
				mImageArr[i*fftLength+j+startIndex]/=weight;
			}
		}

		// norm => 根据fft长度再次归一化
		for(int j=0;j<len;j++){
			mReaArr[i*fftLength+j+startIndex]*=(lenArr[i]/fftLength); 
			mImageArr[i*fftLength+j+startIndex]*=(lenArr[i]/fftLength);
		}

		free(arr);
		free(wArr);
	}
}

// type 0 fftLength/2+1 1 fftLength
static void __cqt_filterBank(int num,float *freBandArr,int samplate,
				int binPerOctave,SpectralFilterBankNormalType normType,WindowType *winType,
				float *factor,float *beta,float *thresh,
				float *lenArr,int fftLength,
				float *mRealFilterBankArr,float *mImageFilterBankArr,int type){
	FFTObj fftObj=NULL;
	int m=0;

	float *mRealTempArr=NULL; // tempKernel
	float *mImageTempArr=NULL;

	float *_mRealArr1=NULL;
	float *_mImageArr1=NULL;

	float _factor=1;
	float _beta=0;
	float _thresh=0.01; // paper 0.0054

	int mLen=0;

	if(!mRealFilterBankArr||!mImageFilterBankArr){
		return;
	}

	if(factor){
		_factor=*factor;
	}

	if(beta){
		_beta=*beta;
	}

	if(thresh){
		_thresh=*thresh;
	}

	m=util_powerTwoBit(fftLength);
	fftObj_new(&fftObj, m);

	_mRealArr1=__vnew(num*fftLength, NULL);
	_mImageArr1=__vnew(num*fftLength, NULL);

	mRealTempArr=__vnew(num*fftLength, NULL);
	mImageTempArr=__vnew(num*fftLength, NULL);

	// 1. tempKernel计算
	__cqt_calTempArr(num,freBandArr,samplate,
				lenArr,fftLength,
				normType,winType,
				mRealTempArr,mImageTempArr);

	// 2. freKernel计算
	for(int i=0;i<num;i++){
		int index=0;

		index=i*fftLength;
		fftObj_fft(fftObj,
				mRealTempArr+index,mImageTempArr+index,
				_mRealArr1+index,_mImageArr1+index);
	}

	// 3.阈值过滤
	if(!type){
		mLen=fftLength/2+1;
	}
	else{
		mLen=fftLength;
	}

	_thresh*=_thresh;
	for(int i=0;i<num;i++){
		for(int j=0;j<mLen;j++){
			float _value1=0;
			float _value2=0;

			_value1=_mRealArr1[i*fftLength+j];
			_value2=_mImageArr1[i*fftLength+j];
			if(_value1*_value1+_value2*_value2>_thresh){
				mRealFilterBankArr[i*mLen+j]=_value1;
				mImageFilterBankArr[i*mLen+j]=_value2;
			}
		}
	}

	free(mRealTempArr);
	free(mImageTempArr);

	free(_mRealArr1);
	free(_mImageArr1);

	fftObj_free(fftObj);
}



























