// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "auditory_weight.h"

void auditory_weightA(float *freArr,int length,float *dBArr){
	float *arr=NULL;

	float min=-80;
	float cArr[]={12200.0*12200.0,20.6*20.6,107.7*107.7,737.9*737.9};

	if(dBArr){
		arr=dBArr;
	}
	else{
		arr=freArr;
	}

	__vpow(freArr, 2, length, arr);
	for(int i=0;i<length;i++){
		float _value=0;

		_value=2.0+20*(log10f(cArr[0])+
				2*log10f(arr[i])-
				log10f(arr[i]+cArr[0])-
				log10f(arr[i]+cArr[1])-
				0.5*log10f(arr[i]+cArr[2])-
				0.5*log10f(arr[i]+cArr[3]));

		arr[i]=(_value<min?min:_value);
	}
}

void auditory_weightB(float *freArr,int length,float *dBArr){
	float *arr=NULL;

	float min=-80;
	float cArr[]={12194.0*12194.0,20.6*20.6,158.5*158.5};

	if(dBArr){
		arr=dBArr;
	}
	else{
		arr=freArr;
	}

	__vpow(freArr, 2, length, arr);
	for(int i=0;i<length;i++){
		float _value=0;

		_value=0.17+20*(log10f(cArr[0])+
				1.5*log10f(arr[i])-
				log10f(arr[i]+cArr[0])-
				log10f(arr[i]+cArr[1])-
				0.5*log10f(arr[i]+cArr[2]));

		arr[i]=(_value<min?min:_value);
	}
}

void auditory_weightC(float *freArr,int length,float *dBArr){
	float *arr=NULL;

	float min=-80;
	float cArr[]={12194.0*12194.0,20.6*20.6};

	if(dBArr){
		arr=dBArr;
	}
	else{
		arr=freArr;
	}

	__vpow(freArr, 2, length, arr);
	for(int i=0;i<length;i++){
		float _value=0;

		_value=0.062+20*(log10f(cArr[0])+
				log10f(arr[i])-
				log10f(arr[i]+cArr[0])-
				log10f(arr[i]+cArr[1]));

		arr[i]=(_value<min?min:_value);
	}
}

void auditory_weightD(float *freArr,int length,float *dBArr){
	float *arr=NULL;

	float min=-80;
	float cArr[]={(8.3046305e-3)*(8.3046305e-3),1018.7*1018.7,1039.6*1039.6,3136.5*3136.5,
					3424.0*3424.0,282.7*282.7,1160.0*1160.0};

	if(dBArr){
		arr=dBArr;
	}
	else{
		arr=freArr;
	}

	__vpow(freArr, 2, length, arr);
	for(int i=0;i<length;i++){
		float _value=0;

		_value=20*(
				0.5*log10f(arr[i])-
				log10f(cArr[0])+
				0.5*(
					log10f((cArr[1]-arr[i])*(cArr[1]-arr[i])+cArr[2]*arr[i])-
					log10f((cArr[3]-arr[i])*(cArr[1]-arr[i])+cArr[4]*arr[i])-
					log10f(cArr[5]+arr[i])-
					log10f(cArr[6]+arr[i])
					));

		arr[i]=(_value<min?min:_value);
	}
}



























