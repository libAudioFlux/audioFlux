// clang -g -c fft_algorithm.c -o fft.o && ar -r libfft.a fft.o

#include <string.h>
#include <math.h>

#ifdef HAVE_ACCELERATE
#include <Accelerate/Accelerate.h>

#elif defined HAVE_FFTW3F
#include <fftw3.h>
#include <pthread.h>

#endif

#include "fft_algorithm.h"

typedef enum{
	FFTExec_None=0,

	FFTExec_FFT,
	FFTExec_IFFT,

	FFTExec_DCT,
	FFTExec_IDCT,

} FFTExecType;

struct OpaqueFFT{
	FFTExecType execType; 

	int radix2Exp; // 2 power
	int fftLength; // length
	int *indexArr; // reverse order cache

	float s0; // dct scale coef
	float s1; 

	int wLength; // w length; fftLength/2

	float *wCosArr; // fft w cache
	float *wSinArr;

	float *wCosArr1; // dct w cache; fftLength
	float *wSinArr1;

	float *realArr1; // fft input r,i cache
	float *imageArr1;

	float *realArr2; // fft output r,i cache
	float *imageArr2;

	int rFlag; // real flag

	#ifdef HAVE_ACCELERATE
	FFTSetup setup;

	DSPSplitComplex inData;
	DSPSplitComplex outData;

	#elif defined HAVE_FFTW3F
	fftwf_plan planRealForward; // real
	fftwf_plan planForward; // complex
	fftwf_plan planInverse; // inverse

	float *inData1;
	fftwf_complex *outData1;

	fftwf_complex *inData2;
	fftwf_complex *outData2;

	fftwf_complex *inData3;
	fftwf_complex *outData3;

	#endif
};

static int *_createIndexArr(int n,int length);

static void _createWArr(float *wCosArr,float *wSinArr,int wLength);
static void _createWArr1(float *wCosArr,float *wSinArr,int wLength);

static void _vmulReal(float *realArr,float *imageArr,float *wCosArr,float *wSinArr,int length);

static void _fftObj_init(FFTObj fftObj);

// flag 0 forward 1 backward
static void _fftObj_fft(FFTObj fftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2,int flag);

int fftObj_new(FFTObj *fftObj,int radix2Exp){
	int status=0;
	FFTObj fft=NULL;

	int length=0;
	int *indexArr=NULL; // reverse order cache

	float s0=0;
	float s1=0;

	int wLength=0;
	float *wCosArr=NULL; // fft w cache
	float *wSinArr=NULL;

	float *wCosArr1=NULL; // dct w cache; length
	float *wSinArr1=NULL;
	
	float *realArr1=NULL;
	float *imageArr1=NULL;

	float *realArr2=NULL;
	float *imageArr2=NULL;

	if(radix2Exp<1||radix2Exp>30){
		status=-100;
		return status;
	}

	fft=*fftObj=(FFTObj )calloc(1, sizeof(struct OpaqueFFT ));

	length=1<<radix2Exp;
	wLength=1<<(radix2Exp-1);

	s0=sqrtf(1.0/length);
	s1=sqrtf(2.0/length);

	realArr1=(float *)calloc(length, sizeof(float ));
	imageArr1=(float *)calloc(length,sizeof(float ));

	realArr2=(float *)calloc(length, sizeof(float ));
	imageArr2=(float *)calloc(length,sizeof(float ));

	// w cache
	wCosArr=(float *)calloc(wLength, sizeof(float ));
	wSinArr=(float *)calloc(wLength, sizeof(float ));
	_createWArr(wCosArr, wSinArr, wLength);

	wCosArr1=(float *)calloc(length, sizeof(float ));
	wSinArr1=(float *)calloc(length, sizeof(float ));
	_createWArr1(wCosArr1, wSinArr1, length);

	// reverse order cache
	indexArr=_createIndexArr(radix2Exp, length);

	fft->radix2Exp=radix2Exp;
	fft->fftLength=length;
	fft->indexArr=indexArr;

	fft->s0=s0;
	fft->s1=s1;

	fft->wLength=wLength;

	fft->wCosArr=wCosArr;
	fft->wSinArr=wSinArr;

	fft->wCosArr1=wCosArr1;
	fft->wSinArr1=wSinArr1;
	
	fft->realArr1=realArr1;
	fft->imageArr1=imageArr1;

	fft->realArr2=realArr2;
	fft->imageArr2=imageArr2;

	_fftObj_init(fft);

	return status;
}

static void _fftObj_init(FFTObj fftObj){
	int radix2Exp=0; // 2 power
	int fftLength=0;

	radix2Exp=fftObj->radix2Exp;
	fftLength=fftObj->fftLength;

	#ifdef HAVE_ACCELERATE
	fftObj->setup=vDSP_create_fftsetup(radix2Exp,FFT_RADIX2);

	#elif defined HAVE_FFTW3F
	fftObj->inData1=(float *)fftwf_malloc(fftLength*sizeof(float ));
	fftObj->outData1=(fftwf_complex *)fftwf_malloc(fftLength*sizeof(fftwf_complex ));

	fftObj->planRealForward=fftwf_plan_dft_r2c_1d(fftLength,
												fftObj->inData1,fftObj->outData1,
												FFTW_ESTIMATE);

	fftObj->inData2=(fftwf_complex *)fftwf_malloc(fftLength*sizeof(fftwf_complex ));
	fftObj->outData2=(fftwf_complex *)fftwf_malloc(fftLength*sizeof(fftwf_complex ));

	fftObj->planForward=fftwf_plan_dft_1d(fftLength,
										fftObj->inData2,fftObj->outData2,
										FFTW_FORWARD,FFTW_ESTIMATE);

	fftObj->inData3=(fftwf_complex *)fftwf_malloc(fftLength*sizeof(fftwf_complex ));
	fftObj->outData3=(fftwf_complex *)fftwf_malloc(fftLength*sizeof(fftwf_complex ));

	fftObj->planInverse=fftwf_plan_dft_1d(fftLength,
										fftObj->inData3,fftObj->outData3,
										FFTW_BACKWARD,FFTW_ESTIMATE);

	#endif

}

#ifdef HAVE_ACCELERATE
static void _fftObj_fft(FFTObj fftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2,int flag){
	int radix2Exp=0;
	int fftLength=0;

	int rFlag=0;

	radix2Exp=fftObj->radix2Exp;
	fftLength=fftObj->fftLength;

	rFlag=fftObj->rFlag;

	fftObj->inData.realp=realArr1;
	fftObj->inData.imagp=imageArr1;

	fftObj->outData.realp=realArr2;
	fftObj->outData.imagp=imageArr2;

	// if(rFlag){ // real to complex fft has problem !!!
	// 	vDSP_fft_zrop(fftObj->setup,
	// 				&fftObj->inData,1,&fftObj->outData,1,
	// 				radix2Exp,FFT_FORWARD);
	// }
	// else{
	// 	vDSP_fft_zop(fftObj->setup,
	// 				&fftObj->inData,1,&fftObj->outData,1,
	// 				radix2Exp,FFT_FORWARD);
	// }

	vDSP_fft_zop(fftObj->setup,
				&fftObj->inData,1,&fftObj->outData,1,
				radix2Exp,flag?FFT_INVERSE:FFT_FORWARD);

}

#elif defined HAVE_FFTW3F
static void _fftObj_fft(FFTObj fftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2,int flag){
	int radix2Exp=0;
	int fftLength=0;

	int rFlag=0;

	radix2Exp=fftObj->radix2Exp;
	fftLength=fftObj->fftLength;

	rFlag=fftObj->rFlag;

	if(rFlag&&!flag){ // real&&forward
		memcpy(fftObj->inData1, realArr1, sizeof(float )*fftLength);

		fftwf_execute(fftObj->planRealForward);

		for(int i=0;i<=fftLength/2;i++){
			realArr2[i]=fftObj->outData1[i][0];
			imageArr2[i]=fftObj->outData1[i][1];
		}

		for(int i=fftLength/2+1,j=fftLength/2-1;i<fftLength;i++,j--){
			realArr2[i]=realArr2[j];
			imageArr2[i]=-imageArr2[j];
		}
	}
	else{ // complex
		if(!flag){ // forward
			for(int i=0;i<fftLength;i++){
				fftObj->inData2[i][0]=realArr1[i];
				fftObj->inData2[i][1]=imageArr1[i];
			}

			fftwf_execute(fftObj->planForward);

			for(int i=0;i<fftLength;i++){
				realArr2[i]=fftObj->outData2[i][0];
				imageArr2[i]=fftObj->outData2[i][1];
			}
		}
		else{ // inverse
			for(int i=0;i<fftLength;i++){
				fftObj->inData3[i][0]=realArr1[i];
				fftObj->inData3[i][1]=imageArr1[i];
			}

			fftwf_execute(fftObj->planInverse);

			for(int i=0;i<fftLength;i++){
				realArr2[i]=fftObj->outData3[i][0];
				imageArr2[i]=fftObj->outData3[i][1];
			}
		}
	}

}

#else
static void _fftObj_fft(FFTObj fftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2,int flag){
	int radix2Exp=0;
	int length=0;

	int *indexArr=NULL; // reverse order cache

	float *wCosArr=NULL; // w cache
	float *wSinArr=NULL;

	int p=0; // 2 power
	int s=0; // step
	int w=0; // w
	int g=0; // group
	int b=0; // butterfly

	float wCos=0;
	float wSin=0;

	float tReal=0;
	float tImage=0;

	radix2Exp=fftObj->radix2Exp;
	length=fftObj->fftLength;

	indexArr=fftObj->indexArr;

	wCosArr=fftObj->wCosArr;
	wSinArr=fftObj->wSinArr;

	// reverse
	for(int i=0;i<length;i++){
		realArr2[i]=realArr1[indexArr[i]];
		imageArr2[i]=imageArr1[indexArr[i]];
	}

	// butterfly algorithm
	for(p=1;p<=radix2Exp;p++){
		s=1<<(p-1);

		for(g=0;g<=s-1;g++){
			w=g*(1<<(radix2Exp-p));
			wCos=wCosArr[w];
			wSin=wSinArr[w];

			for(b=g;b<=length-1;b=b+(1<<p)){
				tReal=realArr2[b+s]*wCos-imageArr2[b+s]*wSin;
				tImage=realArr2[b+s]*wSin+imageArr2[b+s]*wCos;

				realArr2[b+s]=realArr2[b]-tReal;
				imageArr2[b+s]=imageArr2[b]-tImage;

				realArr2[b]=realArr2[b]+tReal;
				imageArr2[b]=imageArr2[b]+tImage;

			}
		}
	}
}

#endif

void fftObj_fft(FFTObj fftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2){
	int length=0;

	float *_realArr1=NULL;
	float *_imageArr1=NULL;

	length=fftObj->fftLength;

	_realArr1=fftObj->realArr1;
	_imageArr1=fftObj->imageArr1;

	if(realArr1&&!imageArr1){
		fftObj->rFlag=1;
	}

	// ??? 
	if(realArr1){
		memcpy(_realArr1, realArr1, sizeof(float )*length);
	}
	else{
		memset(_realArr1, 0, sizeof(float )*length);
	}

	if(imageArr1){
		memcpy(_imageArr1, imageArr1, sizeof(float )*length);
	}
	else{
		memset(_imageArr1, 0, sizeof(float )*length);
	}

	_fftObj_fft(fftObj,_realArr1,_imageArr1,realArr2,imageArr2,0);

	fftObj->execType=FFTExec_FFT;
	fftObj->rFlag=0;
}

void fftObj_ifft(FFTObj fftObj,float *realArr1,float *imageArr1,float *realArr2,float *imageArr2){
	int length=0;

	float *_realArr1=NULL;
	float *_imageArr1=NULL;

	length=fftObj->fftLength;

	_realArr1=fftObj->realArr1;
	_imageArr1=fftObj->imageArr1;

	// if(realArr1&&!imageArr1){
	// 	fftObj->rFlag=1;
	// }

	// ??? 
	if(realArr1){
		memcpy(_realArr1, realArr1, sizeof(float )*length);
	}
	else{
		memset(_realArr1, 0, sizeof(float )*length);
	}

	if(imageArr1){
		memcpy(_imageArr1, imageArr1, sizeof(float )*length);
	}
	else{
		memset(_imageArr1, 0, sizeof(float )*length);
	}

	// #ifdef HAVE_ACCELERATE
	// _fftObj_fft(fftObj,_realArr1,_imageArr1,realArr2,imageArr2,1);

	// #elif defined HAVE_FFTW3F
	// _fftObj_fft(fftObj,_realArr1,_imageArr1,realArr2,imageArr2,1);

	// #else

	// for(int i=0;i<length;i++){
	// 	_imageArr1[i]=-_imageArr1[i];
	// }

	// _fftObj_fft(fftObj,_realArr1,_imageArr1,realArr2,imageArr2,1);

	// for(int i=0;i<length;i++){
	// 	realArr2[i]=realArr2[i]/length;
	// 	imageArr2[i]=-imageArr2[i]/length;
	// }
	
	// #endif

	for(int i=0;i<length;i++){
		_imageArr1[i]=-_imageArr1[i];
	}

	_fftObj_fft(fftObj,_realArr1,_imageArr1,realArr2,imageArr2,0);

	for(int i=0;i<length;i++){
		realArr2[i]=realArr2[i]/length;
		imageArr2[i]=-imageArr2[i]/length;
	}

	fftObj->execType=FFTExec_IFFT;
	fftObj->rFlag=0;
}

void fftObj_dct(FFTObj fftObj,float *dataArr1,float *dataArr2,int isNorm){
	int length=0;

	float *wCosArr1=NULL;
	float *wSinArr1=NULL;

	float *_realArr1=NULL;
	float *_imageArr1=NULL;

	float *_realArr2=NULL;
	float *_imageArr2=NULL;

	length=fftObj->fftLength;

	wCosArr1=fftObj->wCosArr1;
	wSinArr1=fftObj->wSinArr1;

	_realArr1=fftObj->realArr1;
	_imageArr1=fftObj->imageArr1;

	_realArr2=fftObj->realArr2;
	_imageArr2=fftObj->imageArr2;

	fftObj->rFlag=1;

	memset(_imageArr1,0,sizeof(float )*length);
	for(int i=0;i<length/2;i++){
		_realArr1[i]=dataArr1[i*2];
		_realArr1[length-1-i]=dataArr1[i*2+1];
	}

	_fftObj_fft(fftObj,_realArr1,_imageArr1,dataArr2,_imageArr2, 0);
	_vmulReal(dataArr2,_imageArr2,wCosArr1,wSinArr1,length);

	if(isNorm){
		dataArr2[0]=dataArr2[0]*fftObj->s0;
		for(int i=1;i<length;i++){
			dataArr2[i]=dataArr2[i]*fftObj->s1;
		}
	}

	fftObj->execType=FFTExec_DCT;
	fftObj->rFlag=0;
}

void fftObj_idct(FFTObj fftObj,float *dataArr1,float *dataArr2,int isNorm){
	int length=0;

	float *wCosArr1=NULL;
	float *wSinArr1=NULL;

	float *_realArr1=NULL;
	float *_imageArr1=NULL;

	float *_realArr2=NULL;
	float *_imageArr2=NULL;

	length=fftObj->fftLength;

	wCosArr1=fftObj->wCosArr1;
	wSinArr1=fftObj->wSinArr1;

	_realArr1=fftObj->realArr1;
	_imageArr1=fftObj->imageArr1;

	_realArr2=fftObj->realArr2;
	_imageArr2=fftObj->imageArr2;

	if(isNorm){
		dataArr1[0]/=fftObj->s0;
		for(int i=1;i<length;i++){
			dataArr1[i]/=fftObj->s1;
		}
	}

	dataArr1[0]/=2;
	for(int i=0;i<length;i++){
		dataArr1[i]/=fftObj->wLength;
		_realArr1[i]=dataArr1[i]*wCosArr1[i];
		_imageArr1[i]=dataArr1[i]*wSinArr1[i];
	}

	_fftObj_fft(fftObj,_realArr1,_imageArr1,_realArr2,_imageArr2,0);

	for(int i=0;i<length/2;i++){
		dataArr2[i*2]=_realArr2[i];
		dataArr2[i*2+1]=_realArr2[length-1-i];
	}

	fftObj->execType=FFTExec_IDCT;
}

int fftObj_getFFTLength(FFTObj fftObj){

	return fftObj->fftLength;
}

void fftObj_debug(FFTObj fftObj){
	float *realArr2=NULL;
	float *imageArr2=NULL;
	int length=0;

	length=fftObj->fftLength;

	realArr2=fftObj->realArr2;
	imageArr2=fftObj->imageArr2;
	
	if(fftObj->execType!=FFTExec_None){
		if(fftObj->execType==FFTExec_FFT){
			printf("fft ");
		}
		else if(fftObj->execType==FFTExec_IFFT){
			printf("ifft ");
		}
		else if(fftObj->execType==FFTExec_DCT){
			printf("dct ");
		}
		else if(fftObj->execType==FFTExec_IDCT){
			printf("idct ");
		}

		printf("result is: \n");
		for(int i=0;i<length;i++){
			printf("index %d: %f %f\n",i,realArr2[i],imageArr2[i]);
		}
	}
}

void fftObj_free(FFTObj fftObj){
	int *indexArr=NULL; // reverse order cache

	float *wCosArr=NULL; // w cache
	float *wSinArr=NULL;

	float *wCosArr1=NULL; 
	float *wSinArr1=NULL;

	float *realArr1=NULL;
	float *imageArr1=NULL;

	float *realArr2=NULL;
	float *imageArr2=NULL;

	if(!fftObj){
		return;
	}

	indexArr=fftObj->indexArr;

	wCosArr=fftObj->wCosArr;
	wSinArr=fftObj->wSinArr;

	wCosArr1=fftObj->wCosArr1;
	wSinArr1=fftObj->wSinArr1;
	
	realArr1=fftObj->realArr1;
	imageArr1=fftObj->imageArr1;

	realArr2=fftObj->realArr2;
	imageArr2=fftObj->imageArr2;

	free(indexArr);

	free(wCosArr);
	free(wSinArr);

	free(wCosArr1);
	free(wSinArr1);
	
	free(realArr1);
	free(imageArr1);

	free(realArr2);
	free(imageArr2);

	#ifdef HAVE_ACCELERATE
	vDSP_destroy_fftsetup(fftObj->setup);

	#elif defined HAVE_FFTW3F
	fftwf_destroy_plan(fftObj->planRealForward);
	fftwf_destroy_plan(fftObj->planForward);
	fftwf_destroy_plan(fftObj->planInverse);

	fftwf_free(fftObj->inData1);
	fftwf_free(fftObj->outData1);

	fftwf_free(fftObj->inData2);
	fftwf_free(fftObj->outData2);

	fftwf_free(fftObj->inData3);
	fftwf_free(fftObj->outData3);

	#endif

	free(fftObj);
}

static int *_createIndexArr(int n,int length){
	int *indexArr=NULL;
	int *arr;

	int j=0;
	int k=0;

	int _v=0;

	indexArr=(int *)calloc(length,sizeof(int ));
	arr=(int *)calloc(length,sizeof(int ));
	for(int i=0;i<length;i++){
		arr[i]=i;
	}

	for(int i=0;i<length;i++){
		j=0;
		k=0;
		do{
			_v=(i>>(j++))&0x01;
			k=k|_v;
			if(j<n){
				k=k<<1;
			}
		}while(j<n);

		if(k<length){
			indexArr[k]=arr[i];
		}
	}

	free(arr);
	return indexArr;
}

static void _createWArr(float *wCosArr,float *wSinArr,int wLength){
	int length=0;

	length=2*wLength;
	for(int i=0;i<wLength;i++){
		wCosArr[i]=cosf(2*M_PI*i/length);
		wSinArr[i]=-sinf(2*M_PI*i/length);
	}

}

static void _createWArr1(float *wCosArr,float *wSinArr,int wLength){
	int length=0;

	length=2*wLength;
	for(int i=0;i<wLength;i++){
		wCosArr[i]=cosf(M_PI*i/length);
		wSinArr[i]=-sinf(M_PI*i/length);
	}
}

static void _vmulReal(float *realArr,float *imageArr,float *wCosArr,float *wSinArr,int length){

	for(int i=0;i<length;i++){
		realArr[i]=realArr[i]*wCosArr[i]-imageArr[i]*wSinArr[i];
	}
}



