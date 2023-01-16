// clang -g -c 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_complex.h"

#include "flux_window.h"
#include "fft_algorithm.h"
#include "czt_algorithm.h"

struct OpaqueCZT{
	FFTObj fftObj;

	int fftLength; // 2*dataLength

	float lowW; // w=f/fs
	float highW;

	float *nArr; 

	// AW相关
	float *aRealArr; // A^(-n)
	float *aImageArr;

	float *wRealArr; // W^(n^2/2)
	float *wImageArr;

	float *realArr1; // g
	float *imageArr1; 

	float *realArr2; // h
	float *imageArr2; 
	
	// g/h相关 fft
	float *gRealArr;
	float *gImageArr;

	float *hRealArr;
	float *hImageArr;

};

static void _cztObj_initAW(CZTObj cztObj,int fftLength);

static void _cztObj_dealAW(CZTObj cztObj,float lowW,float highW,int fftLength);
static void _cztObj_czt(CZTObj cztObj,float *realArr1,float *imageArr1,float *realArr3,float *imageArr3);

int cztObj_new(CZTObj *cztObj,int radix2Exp){
	int status=0;

	FFTObj fftObj;
	int fftLength=0;

	float *nArr=NULL;

	CZTObj czt=NULL;

	czt=*cztObj=(CZTObj )calloc(1, sizeof(struct OpaqueCZT ));

	status=fftObj_new(&fftObj, radix2Exp+1);
	fftLength=(1<<(radix2Exp+1));

	czt->fftObj=fftObj;
	czt->fftLength=fftLength;

	nArr=__vnew(fftLength, NULL);
	for(int i=-fftLength/2+1,j=0;i<0;i++,j++){
		nArr[j]=i;
	}
	for(int i=0;i<fftLength/2;i++){
		nArr[fftLength/2-1+i]=i;
	}
	czt->nArr=nArr;

	_cztObj_initAW(czt,fftLength);

	return status;
}

void cztObj_czt(CZTObj cztObj,float *realArr1,float *imageArr1,
				float lowW1,float highW1,
				float *realArr3,float *imageArr3){
	
	_cztObj_dealAW(cztObj,lowW1,highW1,cztObj->fftLength);
	_cztObj_czt(cztObj,realArr1,imageArr1,realArr3,imageArr3);

}

void _cztObj_initAW(CZTObj cztObj,int fftLength){

	cztObj->aRealArr=__vnew(fftLength, NULL);
	cztObj->aImageArr=__vnew(fftLength, NULL);

	cztObj->wRealArr=__vnew(fftLength, NULL);
	cztObj->wImageArr=__vnew(fftLength, NULL);

	cztObj->realArr1=__vnew(fftLength, NULL);
	cztObj->imageArr1=__vnew(fftLength, NULL);

	cztObj->realArr2=__vnew(fftLength, NULL);
	cztObj->imageArr2=__vnew(fftLength, NULL);

	cztObj->gRealArr=__vnew(fftLength, NULL);
	cztObj->gImageArr=__vnew(fftLength, NULL);

	cztObj->hRealArr=__vnew(fftLength, NULL);
	cztObj->hImageArr=__vnew(fftLength, NULL);

	_cztObj_dealAW(cztObj,0,1,fftLength);
}

void _cztObj_dealAW(CZTObj cztObj,float lowW1,float highW1,int fftLength){
	float tA=0; // A exp(-j*2*pi*wmin)
	float tW=0; // W exp(-j*2*pi*(wmax-wmin)/N)

	float *nArr=NULL;

	// A/W相关
	float *aRealArr=NULL; // A^(-n)
	float *aImageArr=NULL;

	float *wRealArr=NULL; // W^(n^2/2)
	float *wImageArr=NULL;

	if(lowW1>=highW1||
		lowW1<0||highW1<0||
		highW1>1){
		return;
	}

	if(cztObj->lowW==lowW1&&cztObj->highW==highW1){
		return;
	}

	nArr=cztObj->nArr;
	fftLength=cztObj->fftLength;

	tA=2*M_PI*lowW1;
	tW=-2*M_PI*(highW1-lowW1)/(fftLength/2);

	aRealArr=cztObj->aRealArr;
	aImageArr=cztObj->aImageArr;

	wRealArr=cztObj->wRealArr;
	wImageArr=cztObj->wImageArr;

	for(int i=0;i<fftLength;i++){
		float _n1=0;
		float _n2=0;

		_n1=-nArr[i];
		aRealArr[i]=cosf(_n1*tA);
		aImageArr[i]=sinf(_n1*tA);

		_n2=nArr[i]*nArr[i]/2;
		wRealArr[i]=cosf(_n2*tW);
		wImageArr[i]=sinf(_n2*tW);
	}
}

static void _cztObj_czt(CZTObj cztObj,float *realArr,float *imageArr,float *realArr3,float *imageArr3){
	int fftLength=0;
	FFTObj fftObj=NULL;

	float *aRealArr=NULL; // A^(-n)
	float *aImageArr=NULL;

	float *wRealArr=NULL; // W^(n^2/2)
	float *wImageArr=NULL;

	float *realArr1=NULL; // g
	float *imageArr1=NULL; 

	float *realArr2=NULL; // h
	float *imageArr2=NULL;
	
	// g/h相关 fft
	float *gRealArr=NULL;
	float *gImageArr=NULL;

	float *hRealArr=NULL;
	float *hImageArr=NULL;

	int offset=0;

	fftLength=cztObj->fftLength;
	fftObj=cztObj->fftObj;

	aRealArr=cztObj->aRealArr;
	aImageArr=cztObj->aImageArr;

	wRealArr=cztObj->wRealArr;
	wImageArr=cztObj->wImageArr;

	realArr1=cztObj->realArr1;
	imageArr1=cztObj->imageArr1;

	realArr2=cztObj->realArr2;
	imageArr2=cztObj->imageArr2;

	gRealArr=cztObj->gRealArr;
	gImageArr=cztObj->gImageArr;

	hRealArr=cztObj->hRealArr;
	hImageArr=cztObj->hImageArr;

	// g
	memset(realArr1, 0, sizeof(float )*fftLength);
	memset(imageArr1, 0, sizeof(float )*fftLength);
	offset=fftLength/2-1;
	__vcmul(aRealArr+offset,aImageArr+offset,wRealArr+offset,wImageArr+offset,fftLength/2,realArr1,imageArr1);

	// *data
	if(!realArr||!imageArr){
		for(int i=0;i<fftLength;i++){
			float _v1=0;
			float _v2=0;

			if(realArr){ // real
				_v1=realArr[i];
				realArr1[i]*=_v1;
				imageArr1[i]*=_v1;
			}
			else{ // image
				_v1=imageArr[i];
				_v2=realArr1[i];

				realArr1[i]=-imageArr1[i]*_v1;
				imageArr1[i]=_v2*_v1;
			}
		}
	}
	else{
		__vcmul(realArr1,imageArr1,realArr,imageArr,fftLength,realArr1,imageArr1);
	}
	
	// fft(g)
	fftObj_fft(fftObj,realArr1, imageArr1, gRealArr, gImageArr);

	// fft(h)
	memset(realArr2, 0, sizeof(float )*fftLength);
	memset(imageArr2, 0, sizeof(float )*fftLength);
	memcpy(realArr2, wRealArr, sizeof(float )*(fftLength-1));
	for(int i=0;i<fftLength-1;i++){
		imageArr2[i]=-wImageArr[i];
	}

	fftObj_fft(fftObj,realArr2, imageArr2, hRealArr, hImageArr);

	// result ifft(g*h)*W
	__vcmul(gRealArr, gImageArr, hRealArr, hImageArr,fftLength,gRealArr,gImageArr);
	fftObj_ifft(fftObj, gRealArr,gImageArr,realArr3,imageArr3);
	__vcmul(realArr3+offset, imageArr3+offset, wRealArr+offset, wImageArr+offset,fftLength/2,realArr3,imageArr3);

}

void cztObj_free(CZTObj cztObj){
	FFTObj fftObj=NULL;

	float *aRealArr=NULL; // A^(-n)
	float *aImageArr=NULL;

	float *wRealArr=NULL; // W^(n^2/2)
	float *wImageArr=NULL;

	float *realArr1=NULL; // g
	float *imageArr1=NULL; 

	float *realArr2=NULL; // h
	float *imageArr2=NULL; 
	
	// g/h相关 fft
	float *gRealArr=NULL;
	float *gImageArr=NULL;

	float *hRealArr=NULL;
	float *hImageArr=NULL;

	float *nArr=NULL; 

	if(!cztObj){
		return;
	}

	fftObj=cztObj->fftObj;

	aRealArr=cztObj->aRealArr;
	aImageArr=cztObj->aImageArr;

	wRealArr=cztObj->wRealArr;
	wImageArr=cztObj->wImageArr;

	realArr1=cztObj->realArr1;
	imageArr1=cztObj->imageArr1;

	realArr2=cztObj->realArr2;
	imageArr2=cztObj->imageArr2;

	gRealArr=cztObj->gRealArr;
	gImageArr=cztObj->gImageArr;

	hRealArr=cztObj->hRealArr;
	hImageArr=cztObj->hImageArr;

	nArr=cztObj->nArr;

	fftObj_free(fftObj);

	free(nArr);

	free(aRealArr);
	free(aImageArr);

	free(wRealArr);
	free(wImageArr);

	free(realArr1);
	free(imageArr1);

	free(realArr2);
	free(imageArr2);

	free(gRealArr);
	free(gImageArr);

	free(hRealArr);
	free(hImageArr);

	free(cztObj);
}






