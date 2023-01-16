// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_complex.h"

#include "filterDesign_freqz.h"

// isWhole 0默认计算一半
void filterDesign_freqzBA(float *bArr,float *aArr,int length,
						int fftLength,int samplate,int isWhole,
						float *kArr,
						float *realArr,float *imageArr,
						float *wArr){
	float *rArr=NULL;
	float *iArr=NULL;

	float end=0;
	float *_kArr=NULL;

	int len=0;
	
	end=2*M_PI;
	if(kArr){
		_kArr=kArr;
	}
	else{
		_kArr=__vlinspace(0, end-end/fftLength, fftLength, 0);
	}

	if(!isWhole){ // 只计算一半
		len=fftLength/2+1;
	}
	else{
		len=fftLength;
	}

	__vcnew(len, NULL, &rArr, &iArr);

	filterDesign_calFreResponse(_kArr, len, bArr, length, realArr, imageArr);
	filterDesign_calFreResponse(_kArr, len, aArr, length, rArr, iArr);
	__vcdiv(realArr, imageArr, rArr, iArr, len, NULL, NULL);

	if(wArr){
		for(int i=0;i<len;i++){
			wArr[i]=_kArr[i]*samplate/(2*M_PI);
		}
	}

	free(rArr);
	free(iArr);

	if(!kArr){
		free(_kArr);
	}
}

/***
	SOS => second order of section
	传输函数二阶相乘 
	满足二阶mLength=6
****/
void filterDesign_freqzSOS(float *mArr,int nLength,
						int fftLength,int samplate,int isWhole,
						float *kArr,
						float *realArr,float *imageArr,
						float *wArr){
	float *rArr=NULL;
	float *iArr=NULL;

	float end=0;
	float *_kArr=NULL;

	int mLength=0;
	int len=0;

	mLength=6;

	end=2*M_PI;
	if(kArr){
		_kArr=kArr;
	}
	else{
		_kArr=__vlinspace(0, end-end/fftLength, fftLength, 0);
	}

	if(!isWhole){ // 只计算一半
		len=fftLength/2+1;
	}
	else{
		len=fftLength;
	}

	filterDesign_freqzBA(mArr, mArr+mLength/2, mLength/2,
						 fftLength, samplate,isWhole, 
						 _kArr,
						 realArr, imageArr,
						 wArr);

	__vcnew(len, NULL, &rArr, &iArr);
	for(int i=1;i<nLength;i++){	
		filterDesign_freqzBA(mArr+i*mLength, mArr+(i*mLength+mLength/2), mLength/2,
						 fftLength, samplate,isWhole, 
						 _kArr,
						 rArr, iArr,
						 NULL);

		__vcmul(realArr, imageArr, rArr, iArr, len, NULL, NULL);
	}

	free(rArr);
	free(iArr);

	if(!kArr){
		free(_kArr);
	}
}

/***
	b[0]+b[1]*e^(-jw)+...+b[M]*e^(-j(M-1)w);
	wArr1 0~2*pi/0~pi区间角频率刻度 fftLength||fftLength/2
	vArr2 b/a系数
****/
void filterDesign_calFreResponse(float *wArr1,int length1,
							float *vArr2,int length2,
							float *realArr3,float *imageArr3){
	
	for(int i=0;i<length1;i++){
		realArr3[i]=0;
		imageArr3[i]=0;
		for(int j=0;j<length2;j++){
			realArr3[i]+=cosf(-wArr1[i]*j)*vArr2[j];
			imageArr3[i]+=sinf(-wArr1[i]*j)*vArr2[j];
		}
	}
}


















