// 

#include <string.h>
#include <math.h>

#include "flux_vector.h"
#include "flux_vectorOp.h"

// element math相关
// abs/ng/floor/ceil/round
void __vabs(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=fabsf(vArr1[i]);
	}
}

void __vng(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=-vArr1[i];
	}
}

void __vfloor(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=floorf(vArr1[i]);
	}
}

void __vceil(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=ceilf(vArr1[i]);
	}
}

void __vround(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=roundf(vArr1[i]);
	}
}


// cos/sin/tan/acos/asin/atan
void __vcos(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=cosf(vArr1[i]);
	}
}

void __vsin(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=sinf(vArr1[i]);
	}
}

void __vtan(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=tanf(vArr1[i]);
	}	
}	

void __vacos(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=acosf(vArr1[i]);
	}
}

void __vasin(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=asinf(vArr1[i]);
	}
}

void __vatan(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=atanf(vArr1[i]);
	}
}

// exp/exp2/pow/sqrt
void __vexp(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=expf(vArr1[i]);
	}
}

void __vexp2(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=exp2f(vArr1[i]);
	}
}

void __vpow(float *vArr1,float exp,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=powf(vArr1[i],exp);
	}
}

void __vsqrt(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=sqrtf(vArr1[i]);
	}
}


// log/log10/log2
void __vlog(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=logf(vArr1[i]);
	}
}

void __vlog10(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=log10f(vArr1[i]);
	}
}

void __vlog2(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=log2f(vArr1[i]);
	}
}

// log compress
void __vlog_compress(float *vArr1,float *gamma,float *base,int length,float *vArr2){
	float *arr=NULL;

	float _gamma=1.0;
	float _base=1.0;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	if(gamma){
		if(*gamma>0){
			_gamma=*gamma;
		}
	}

	if(base){
		if(*base>0){
			_base=*base;
		}
	}

	for(int i=0;i<length;i++){
		arr[i]=logf(vArr1[i]*_gamma+_base);
	}
}

void __vlog10_compress(float *vArr1,float *gamma,float *base,int length,float *vArr2){
	float *arr=NULL;

	float _gamma=1.0;
	float _base=1.0;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	if(gamma){
		if(*gamma>0){
			_gamma=*gamma;
		}
	}

	if(base){
		if(*base>0){
			_base=*base;
		}
	}

	for(int i=0;i<length;i++){
		arr[i]=log10f(vArr1[i]*_gamma+_base);
	}
}

void __vlog2_compress(float *vArr1,float *gamma,float *base,int length,float *vArr2){
	float *arr=NULL;

	float _gamma=1.0;
	float _base=1.0;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	if(gamma){
		if(*gamma>0){
			_gamma=*gamma;
		}
	}

	if(base){
		if(*base>0){
			_base=*base;
		}
	}

	for(int i=0;i<length;i++){
		arr[i]=log2f(vArr1[i]*_gamma+_base);
	}
}

// sinc=sin(pi*x)/x
void __vsinc(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		float _value=0;

		_value=vArr1[i]*M_PI;
		if(fabsf(_value)<1e-9){
			arr[i]=1;
		}
		else{
			arr[i]=sinf(_value)/_value;
		}
	}

}

void __vsinc_low(float *vArr1,int length,float cut,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	if(cut<=0||cut>=1){
		return;
	}

	for(int i=0;i<length;i++){
		arr[i]*=cut;
	}

	__vsinc(arr,length,arr);

	for(int i=0;i<length;i++){
		arr[i]*=cut;
	}
}

// sinc(1)-sinc(cut)
void __vsinc_high(float *vArr1,int length,float cut,float *vArr2){
	float *arr=NULL;
	float *_arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	if(cut<=0||cut>=1){
		return;
	}

	if(length%2==0){ // high coeff must odd
		printf("high coeff must odd!!!\n");
		return;
	}

	_arr=__vnew(length, NULL);

	__vsinc(vArr1,length,_arr); // sinc(1)
	__vsinc_low(vArr1,length,cut,vArr2); // sinc(cut)

	for(int i=0;i<length;i++){
		arr[i]=_arr[i]-vArr2[i];
	}

	free(_arr);
}

// sinc(cut2)-sinc(cut1)
void __vsinc_bandpass(float *vArr1,int length,float cut1,float cut2,float *vArr2){
	float *arr=NULL;
	float *_arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	if(cut1<=0||cut1>=1){
		return;
	}

	if(cut2<=0||cut2>=1){
		return;
	}

	if(cut2<=cut1){
		return;
	}

	_arr=__vnew(length, NULL);
	for(int i=0;i<length;i++){
		_arr[i]=vArr1[i]; // *cut2
	}
	__vsinc_low(_arr,length,cut2,_arr); // sinc(cut2)

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]; // *cut1
	}
	__vsinc_low(arr,length,cut1,arr); // sinc(cut1)

	for(int i=0;i<length;i++){
		arr[i]=_arr[i]-arr[i];
	}

	free(_arr);
}

// sinc(1)-sinc(cut2)+sinc(cut1)
void __vsinc_stop(float *vArr1,int length,float cut1,float cut2,float *vArr2){
	float *arr=NULL;

	float *_arr1=NULL;
	float *_arr2=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	if(cut1<=0||cut1>=1){
		return;
	}

	if(cut2<=0||cut2>=1){
		return;
	}

	if(cut2<=cut1){
		return;
	}

	if(length%2==0){ // high coeff must odd
		printf("stop coeff must odd!!!\n");
		return;
	}

	_arr1=__vnew(length, NULL);
	_arr2=__vnew(length, NULL);
	for(int i=0;i<length;i++){
		_arr1[i]=vArr1[i]; // *cut2
	}
	__vsinc_low(_arr1,length,cut2,_arr1); // sinc(cut2)

	for(int i=0;i<length;i++){
		_arr2[i]=vArr1[i]; // *cut1
	}
	__vsinc_low(_arr2,length,cut1,_arr2); // sinc(cut1)

	__vsinc(vArr1,length,arr); // sinc(1)
	for(int i=0;i<length;i++){
		arr[i]=arr[i]-_arr1[i]+_arr2[i];
	}

	free(_arr1);
	free(_arr2);
}

// 1/2阶 
void __vgradient(float *vArr1,int length,int edgeOrder,float *vArr2){

	if(length<2){
		return;
	}

	for(int i=0;i<length;i++){
		if(i==0){ // start f(k+1)-f(k)
			vArr2[i]=vArr1[i+1]-vArr1[i];
		}
		else if(i==length-1){ // f(k)-f(k-1)
			vArr2[i]=vArr1[i]-vArr1[i-1];
		}
		else{ // (f(k+1)-f(k-1))/2
			vArr2[i]=(vArr1[i+1]-vArr1[i-1])/2;
		}
	}

	if(edgeOrder>1){ // 2阶
		vArr2[0]-=(vArr2[1]-vArr2[0]);
		vArr2[length-1]+=(vArr2[length-1]-vArr2[length-2]);
	}
}

// interpolation
void __vinterp_linear(float *xArr1,float *yArr1,int length1,float *xArr2,int length2,float *yArr2){
	int index1=0;

	float x1=0;
	float y1=0;

	float x2=0;
	float y2=0;

	float value=0;

	for(int i=0;i<length2;i++){
		while(index1<length1-1&&xArr2[i]>xArr1[index1+1]){
			index1++;
		}

		if(index1<length1-1){
			x1=xArr1[index1];
			y1=yArr1[index1];

			x2=xArr1[index1+1];
			y2=yArr1[index1+1];

			value=y1+(xArr2[i]-x1)*(y2-y1)/(x2-x1);
			yArr2[i]=value;
		}
		else{
			yArr2[i]=yArr1[length1-1];
		}
	}
}

// pad type 0 'constant' 1 'refect' 2 'wrap'
void __vpad(float *vArr1,int length,
			int headLen,int tailLen,int type,
			float *value1,float *value2,
			float *vArr2){

}

void __mpad(float *mArr1,int nlength,int mLength,
			int headLen,int tailLen,int type,
			float *value1,float *value2,
			int axis,
			float *mArr2){

}

// constant
void __vpad_center1(float *vArr1,int vLength,
				int leftLength,int rightLength,
				float value1,float value2){
	for(int i=0;i<leftLength;i++){
		vArr1[i]=value1;
	}

	for(int i=leftLength+vLength;i<leftLength+rightLength+vLength;i++){
		vArr1[i]=value2;
	}
}

void __vpad_left1(float *vArr1,int vLength,int fLength,int value){
	for(int i=0;i<fLength;i++){
		vArr1[i]=value;
	}
}

void __vpad_right1(float *vArr1,int vLength,int fLength,int value){
	for(int i=vLength;i<fLength+vLength;i++){
		vArr1[i]=value;
	}
}

// reflect
void __vpad_center2(float *vArr1,int vLength,int leftLength,int rightLength){
	int startIndex=0;

	int leftIndex=0;
	int rightIndex=0;

	int flag=0; // 0 升 1降

	if(vLength<2){
		return;
	}

	startIndex=leftLength;
	leftIndex=1;
	for(int i=leftLength-1;i>=0;i--){
		vArr1[i]=vArr1[startIndex+leftIndex];

		if(!flag){ // 升
			if(leftIndex==vLength-1){
				flag=!flag;
				leftIndex--;
				continue;
			}
			else{
				leftIndex++;
			}
		}
		else{
			if(leftIndex==0){
				flag=!flag;
				leftIndex++;
				continue;
			}
			leftIndex--;
		}

		if(leftIndex==0||leftIndex==vLength-1){ // start||end
			flag=!flag;
		}
	}

	flag=1; // 降
	rightIndex=vLength-2;
	for(int i=leftLength+vLength;i<leftLength+rightLength+vLength;i++){
		vArr1[i]=vArr1[startIndex+rightIndex];

		if(!flag){ // 升
			if(rightIndex==vLength-1){
				flag=!flag;
				rightIndex--;
				continue;
			}
			else{
				rightIndex++;
			}
		}
		else{
			if(rightIndex==0){
				flag=!flag;
				rightIndex++;
				continue;
			}
			rightIndex--;
		}

		if(rightIndex==0||rightIndex==vLength-1){ // start||end
			flag=!flag;
		}
	}
}

void __vpad_left2(float *vArr1,int vLength,int leftLength){
	
	__vpad_center2(vArr1,vLength,leftLength,0);
}

void __vpad_right2(float *vArr1,int vLength,int rightLength){
	
	__vpad_center2(vArr1,vLength,0,rightLength);
}

// wrap
void __vpad_center3(float *vArr1,int vLength,int leftLength,int rightLength){
	int startIndex=0;

	int leftIndex=0;
	int rightIndex=0;

	if(vLength<2){
		return;
	}

	startIndex=leftLength;
	leftIndex=vLength-1;
	for(int i=leftLength-1;i>=0;i--){
		vArr1[i]=vArr1[startIndex+leftIndex];

		if(leftIndex==0){
			leftIndex=vLength-1;
		}
		else{
			leftIndex--;
		}
	}

	rightIndex=0;
	for(int i=leftLength+vLength;i<leftLength+rightLength+vLength;i++){
		vArr1[i]=vArr1[startIndex+rightIndex];

		if(rightIndex==vLength-1){
			rightIndex=0;
		}
		else{
			rightIndex++;
		}
	}
}

void __vpad_left3(float *vArr1,int vLength,int leftLength){
	
	__vpad_center3(vArr1,vLength,leftLength,0);
}

void __vpad_right3(float *vArr1,int vLength,int rightLength){

	__vpad_center3(vArr1,vLength,0,rightLength);
}

// symmetric 和reflect接近 min/max
void __vpad_center4(float *vArr1,int vLength,int leftLength,int rightLength){
	int startIndex=0;

	int leftIndex=0;
	int rightIndex=0;

	int flag=0; // 0 升 1降

	if(vLength<2){
		return;
	}

	startIndex=leftLength;
	leftIndex=0;
	for(int i=leftLength-1;i>=0;i--){
		vArr1[i]=vArr1[startIndex+leftIndex];

		if(!flag){ // 升
			if(leftIndex==vLength-1){
				flag=!flag;
				continue;
			}
			else{
				leftIndex++;
			}
		}
		else{
			if(leftIndex==0){
				flag=!flag;
				continue;
			}
			leftIndex--;
		}
	}

	flag=1; // 降
	rightIndex=vLength-1;
	for(int i=leftLength+vLength;i<leftLength+rightLength+vLength;i++){
		vArr1[i]=vArr1[startIndex+rightIndex];

		if(!flag){ // 升
			if(rightIndex==vLength-1){
				flag=!flag;
				continue;
			}
			else{
				rightIndex++;
			}
		}
		else{
			if(rightIndex==0){
				flag=!flag;
				continue;
			}
			rightIndex--;
		}
	}
}

void __vpad_left4(float *vArr1,int vLength,int leftLength){

	__vpad_center4(vArr1,vLength,leftLength,0);
}

void __vpad_right4(float *vArr1,int vLength,int rightLength){

	__vpad_center4(vArr1,vLength,0,rightLength);
}

// edge
void __vpad_center5(float *vArr1,int vLength,int leftLength,int rightLength){
	for(int i=0;i<leftLength;i++){
		vArr1[i]=vArr1[leftLength];
	}

	for(int i=leftLength+vLength;i<leftLength+rightLength+vLength;i++){
		vArr1[i]=vArr1[leftLength+vLength-1];
	}
}

void __vpad_left5(float *vArr1,int vLength,int leftLength){
	for(int i=0;i<leftLength;i++){
		vArr1[i]=vArr1[leftLength];
	}
}

void __vpad_right5(float *vArr1,int vLength,int rightLength){
	for(int i=vLength;i<vLength+rightLength;i++){
		vArr1[i]=vArr1[vLength-1];
	}
}









