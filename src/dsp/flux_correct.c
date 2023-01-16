// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "flux_correct.h"

void correct_rect(float cur,float left,float right,float *det,float *value){
	float _det=0;
	float _value=0;

	double eps=1e-10;

	float y1=0;
	float y2=0;

	float v1=0;
	float v2=0;

	float n=0;
	float s=0;

	float c1=0;
	float c2=0;

	if(right>=left){
		y1=cur;
		y2=right;
	}
	else{
		y1=left;
		y2=cur;
	}

	if(y2<eps){
		y2=eps;
	}

	// 1. det计算
	v1=y1/y2;
	v2=1+v1;
	if(v2<eps){
		v2=eps;
	}

	_det=1/v2;
	if(y1<y2){
		_det-=1;
	}

	// 2. value计算
	if(_det>=0){
		n=floorf(_det);
	}
	else{
		n=ceilf(_det);
	}
	s=_det-n;
	
	if(fabs(s)<1e-8){
		s=1e-8;
	}

	c1=n+s;
	c2=M_PI*c1/sinf(M_PI*c1);

	_value=cur*c2;

	if(det){
		*det=_det;
	}

	if(value){
		*value=_value;
	}
}

void correct_hann(float cur,float left,float right,float *det,float *value){
	float _det=0;
	float _value=0;

	double eps=1e-10;

	float y1=0;
	float y2=0;

	float v1=0;
	float v2=0;

	float n=0;
	float s=0;

	float c1=0;
	float c2=0;

	if(right>=left){
		y1=cur;
		y2=right;
	}
	else{
		y1=left;
		y2=cur;
	}

	if(y2<eps){
		y2=eps;
	}

	// 1. det计算
	v1=y1/y2;
	v2=1+v1;
	if(v2<eps){
		v2=eps;
	}

	_det=(2-v1)/v2;
	if(y1<y2){
		_det-=1;
	}

	// 2. value计算
	if(_det>=0){
		n=floorf(_det);
	}
	else{
		n=ceilf(_det);
	}
	s=_det-n;
	
	if(fabs(s)<1e-8){
		s=1e-8;
	}

	c1=n+s;
	c2=M_PI*c1/sinf(M_PI*c1);

	_value=cur*c2*(1-c1*c1)*2;

	if(det){
		*det=_det;
	}

	if(value){
		*value=_value;
	}
}

void correct_hamm(float cur,float left,float right,float *det,float *value){
	float _det=0;
	float _value=0;

	double eps=1e-10;

	float y1=0;
	float y2=0;

	float v1=0;
	float v2=0;

	float n=0;
	float s=0;

	float c1=0;
	float c2=0;

	if(right>=left){
		y1=cur;
		y2=right;
	}
	else{
		y1=left;
		y2=cur;
	}

	if(y2<eps){
		y2=eps;
	}

	// 1. det计算
	c1=-27.0/4;
	v1=y1/y2;

	_det=-(2-v1)/(1+v1);
	for(int i=0;i<8;i++){
		v2=(_det*_det+c1)/((_det+1)*(_det+1)+c1);
		_det=(v1-2*v2)/(v1+v2);
	}
	_det=-_det;

	if(y1<y2){
		_det-=1;
	}

	// 2. value计算
	if(_det>=0){
		n=floorf(_det);
	}
	else{
		n=ceilf(_det);
	}
	s=_det-n;
	
	if(fabs(s)<1e-8){
		s=1e-8;
	}

	c1=n+s;
	c2=M_PI*c1/sinf(M_PI*c1);

	_value=cur*c2*(1-c1*c1)/(0.54-0.08*c1*c1);

	if(det){
		*det=_det;
	}

	if(value){
		*value=_value;
	}
}

float correct_getRectRecover(){

	return 1;
}

float correct_getHannRecover(){

	return 1/0.5;
}

float correct_getHammRecover(){

	return 1/0.54;
}










