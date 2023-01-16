// clang -g 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "flux_window.h"
#include "filterDesign_fir.h"

static void __scaleFilter(float *bArr1,int first,float gain,int length,float *bArr2);

/***
	fir 窗函数法设计
	1. 根据滤波器类型计算理想sinc采样
	2. 加窗
	3. 归一化处理
	wc=(ws+wp)/2 wc1<wc2 =>bandpass/stop
****/
float *filterDesign_fir2(int order,
					float *wcArr,
					FilterBandType bandType,
					float *winArr,
					int *isNoScale){
	float *bArr=NULL;
	int isScale=1;

	int first=0;
	float gain=0;


	if(bandType==FilterBand_HighPass||
		bandType==FilterBand_BandStop){ // order must even
		if(order%2!=0){
			printf("high or stop order must even!!!\n");
			return NULL;
		}
	}

	bArr=__vlinspace(-order/2, order/2, order+1, 0);
	if(bandType==FilterBand_LowPass){ // lowpass
		__vsinc_low(bArr,order+1,wcArr[0],bArr);
	}
	else if(bandType==FilterBand_HighPass){ // highpass
		__vsinc_high(bArr,order+1,wcArr[0],bArr);
	}
	else if(bandType==FilterBand_BandPass){ 
		__vsinc_bandpass(bArr,order+1,wcArr[0],wcArr[1],bArr);
	}
	else if(bandType==FilterBand_BandStop){
		__vsinc_stop(bArr,order+1,wcArr[0],wcArr[1],bArr);
	}

	__vmul(bArr, winArr, order+1, bArr);

	__vdebug(bArr, order+1, 1);
	printf("\n\n");

	if(isNoScale){
		if(*isNoScale!=0){
			isScale=0;
		}
	}

	if(isScale){
		if(bandType==FilterBand_LowPass||
			bandType==FilterBand_BandStop){
			first=1;
		}

		if(bandType==FilterBand_HighPass){
			gain=1.0;
		}
		else if(bandType==FilterBand_BandPass){
			gain=(wcArr[0]+wcArr[1])/2;
		}

		__scaleFilter(bArr,first,gain,order+1,bArr);
	}
	
	return bArr;
}

float *filterDesign_fir1(int order,
					float *wcArr,
					FilterBandType bandType,
					WindowType *winType,float *value,
					int *isNoScale){
	float *bArr=NULL;
	float *winArr=NULL;

	WindowType type=Window_Hamm;

	if(bandType==FilterBand_HighPass||
		bandType==FilterBand_BandStop){ // order must even
		if(order%2!=0){
			printf("high or stop order must even!!!\n");
			return NULL;
		}
	}

	if(winType){
		type=*winType;
	}

	if(type==Window_Rect){
		float _value=1;

		winArr=__vnew(order+1, &_value);
	}
	else if(type==Window_Hann){
		winArr=window_createHann(order+1,0);
	}
	else if(type==Window_Hamm){
		winArr=window_createHamm(order+1,0);
	}
	else if(type==Window_Blackman){
		winArr=window_createBlackman(order+1,0);
	}
	else if(type==Window_Kaiser){
		winArr=window_createKaiser(order+1,0,value);
	}
	else if(type==Window_Bartlett){
		winArr=window_createBartlett(order+1,0);
	}
	else if(type==Window_Triang){
		winArr=window_createTriang(order+1,0);
	}
	else if(type==Window_Flattop){
		winArr=window_createFlattop(order+1,0);
	}
	else if(type==Window_Gauss){
		winArr=window_createGauss(order+1,0,value);
	}
	else if(type==Window_Blackman_Harris){
		winArr=window_createBlackmanHarris(order+1,0);
	}
	else if(type==Window_Blackman_Nuttall){
		winArr=window_createBlackmanNuttall(order+1,0);
	}
	else if(type==Window_Bartlett_Hann){
		winArr=window_createBartlettHann(order+1,0);
	}
	else if(type==Window_Bohman){
		winArr=window_createBohman(order+1,0);
	}
	else if(type==Window_Tukey){
		winArr=window_createTukey(order+1,0,value);
	}

	bArr=filterDesign_fir2(order,wcArr,bandType,winArr,isNoScale);

	free(winArr);
	return bArr;
}

static void __scaleFilter(float *bArr1,int first,float gain,int length,float *bArr2){
	float *arr=NULL;

	float rValue=0;
	float iValue=0;

	float value=0;

	if(bArr2){
		arr=bArr2;
	}
	else{
		arr=bArr1;
	}

	if(first){ // lowpass/stop
		value=__vsum(bArr1, length);
		for(int i=0;i<length;i++){
			arr[i]=bArr1[i]/value;
		}
	}
	else{ // high/bandpass
		for(int i=0;i<length;i++){
			rValue+=cosf(2*M_PI*i*(gain/2))*bArr1[i];
			iValue+=-sinf(2*M_PI*i*(gain/2))*bArr1[i];
		}

		value=sqrtf(rValue*rValue+iValue*iValue);
		for(int i=0;i<length;i++){
			arr[i]=bArr1[i]/value;
		}
	}
}

// order must odd; first derivative
float *filterDesign_smooth1(int order){
	float *arr1=NULL;
	
	float v1=0;
	int m=0;

	if(!(order&1)){
		printf("order must odd !!!\n");
		return NULL;
	}

	arr1=__vnew(order+1, NULL);

	m=floorf(order/2);
	for(int i=1;i<=m;i++){
		v1+=i*i;
	}

	for(int i=m,j=0;i>=-m;i--,j++){
		arr1[j]=i/v1;
	}

	return arr1;
}

float *filterDesign_mean(int order){
	float *arr1=NULL;
	float _value=0;

	_value=1.0/order;
	arr1=__vnew(order, &_value);

	return arr1;
}

void filterDesign_filter(float *bArr,float *aArr,float *xArr,
						int bLength,int aLength,int xLength,
						float *yArr){

	yArr[0]=bArr[0]*xArr[0];
	for(int i=1;i<xLength;i++){

		for(int j=0;j<bLength;j++){
			if(i>=j){
				yArr[i]=yArr[i]+bArr[j]*xArr[i-j];
			}
		}

		for(int k=0;k<aLength-1;k++){
			if(i>k){
				yArr[i]=yArr[i]-aArr[k+1]*yArr[i-k-1];
			}
		}
	}
}

// filtfilt(b,a,x)
void filterDesign_filtfilt(float *bArr,float *aArr,float *xArr,
						int bLength,int aLength,int xLength,
						float *yArr){
	
	
}









