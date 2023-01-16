// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"

#include "flux_window.h"

static float *__calHalfHann(int halfLen,int length);
static float *__calHalfHamm(int halfLen,int length);
static float *__calHalfBlackman(int halfLen,int length);
static float *__calHalfFlattop(int halfLen,int length);

static float *__calHalfBartlett(int halfLen,int length);
static float *__calHalfTriang(int halfLen,int length);

static float __besselZeroOne(float value);
// alpha 默认5
static float *__calHalfKaiser(int halfLen,float alpha,int length);
// alpha 默认2.5
static float *__calHalfGauss(int halfLen,float alpha,int length);

static float *__calHalfBlackmanHarris(int halfLen,int length);
static float *__calHalfBlackmanNuttall(int halfLen,int length);
static float *__calHalfBartlettHann(int halfLen,int length);

static float *__calHalfBohman(int halfLen,int length);

static void __calKaiserOrder(float w1,float w2,float atten,int *order,float *beta);

static float *__window_createHann(int length,int flag);
static float *__window_createHamm(int length,int flag);

static float *__window_createBlackman(int length,int flag);
// beta 5.0
static float *__window_createKaiser(int length,float *beta);

// only symmetric
static float *__window_createBartlett(int length,int flag);
static float *__window_createTriang(int length,int flag);

static float *__window_createFlattop(int length,int flag);
// alpha 2.5
static float *__window_createGauss(int length,float *alpha);

// 4-term blackman harris~nuttall
static float *__window_createBlackmanHarris(int length,int flag);
static float *__window_createBlackmanNuttall(int length,int flag);
// like Bartlett,hann,hamm; only symmetric
static float *__window_createBartlettHann(int length,int flag);

// only symmetric
static float *__window_createBohman(int length,int flag);

// alpha>=0&&<=1 0.5 rect~hann ;only symmetric
static float *__window_createTukey(int length,float *alpha);

// order相关 
// w1 ws/wp w2 ws/wp atten dB
void window_calKaiserOrder(float w1,float w2,float atten,int *order,float *beta);

// 窗函数相关 flag 0 symmetric对称 1 periodic周期 针对fft length extend 1 sample
float *window_createHann(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createHann(flag?length+1:length,0);
		}
	}

	return arr;
}

float *window_createHamm(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createHamm(flag?length+1:length,0);
		}
	}

	return arr;
}


float *window_createBlackman(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createBlackman(flag?length+1:length,0);
		}
	}

	return arr;
}

// beta 5.0
float *window_createKaiser(int length,int flag,float *beta){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createKaiser(flag?length+1:length,beta);
		}
	}

	return arr;
}

// only symmetric
float *window_createBartlett(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createBartlett(flag?length+1:length,0);
		}
	}

	return arr;
}

float *window_createTriang(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createTriang(flag?length+1:length,0);
		}
	}

	return arr;
}

float *window_createFlattop(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createFlattop(flag?length+1:length,0);
		}
	}

	return arr;
}

// alpha 2.5
float *window_createGauss(int length,int flag,float *alpha){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createGauss(flag?length+1:length,alpha);
		}
	}

	return arr;
}


// 4-term blackman harris~nuttall
float *window_createBlackmanHarris(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createBlackmanHarris(flag?length+1:length,0);
		}
	}

	return arr;
}

float *window_createBlackmanNuttall(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createBlackmanNuttall(flag?length+1:length,0);
		}
	}

	return arr;
}

// like Bartlett,hann,hamm; only symmetric
float *window_createBartlettHann(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createBartlettHann(flag?length+1:length,0);
		}
	}

	return arr;
}

// only symmetric
float *window_createBohman(int length,int flag){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createBohman(flag?length+1:length,0);
		}
	}

	return arr;
}

// alpha>=0&&<=1 0.5 rect~hann ;only symmetric
float *window_createTukey(int length,int flag,float *alpha){
	float *arr=NULL;

	if(length>0){
		if(length==1){
			arr=__vnew(1, NULL);
			arr[0]=1;
		}
		else{
			arr=__window_createTukey(flag?length+1:length,alpha);
		}
	}

	return arr;
}

static float *__window_createHann(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<1){
		return NULL;
	}
	
	if(!(length&1)){ // even
		halfLen=length/2+flag; // 针对periodic +1
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfHann(halfLen,length-1+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createHamm(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<1){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag; // 针对periodic +1
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfHamm(halfLen,length-1+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createBlackman(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<1){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag; // 针对periodic +1
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfBlackman(halfLen,length-1+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createBlackmanHarris(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<1){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfBlackmanHarris(halfLen,length-1+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createBlackmanNuttall(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<1){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfBlackmanNuttall(halfLen,length-1+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createKaiser(int length,float *alpha){
	float *wArr=NULL;
	int halfLen=0;

	float a=5.0;

	if(length<2){
		return NULL;
	}

	if(alpha){
		if(*alpha>0){
			a=*alpha;
		}
	}

	if(!(length&1)){ // even
		halfLen=length/2;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfKaiser(halfLen,a,length-1);
	for(int i=length-1,j=0;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

// bartlett/triang 三角窗相关 length=N+1 !!!
static float *__window_createBartlett(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<2){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfBartlett(halfLen,length-1+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createBartlettHann(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<2){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfBartlettHann(halfLen,length-1+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createTriang(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<1){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfTriang(halfLen,length+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createFlattop(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<1){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfFlattop(halfLen,length-1+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createGauss(int length,float *alpha){
	float *wArr=NULL;
	int halfLen=0;

	float a=2.5;

	if(length<2){
		return NULL;
	}

	if(alpha){
		if(*alpha>0){
			a=*alpha;
		}
	}

	if(!(length&1)){ // even
		halfLen=length/2+1;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfGauss(halfLen,a,length);
	for(int i=length-1,j=0;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createBohman(int length,int flag){
	float *wArr=NULL;
	int halfLen=0;

	if(length<1){
		return NULL;
	}

	if(!(length&1)){ // even
		halfLen=length/2+flag;
	}
	else{ // odd
		halfLen=(length+1)/2;
	}

	wArr=__calHalfBohman(halfLen,length+flag);
	for(int i=length-1,j=flag;i>=halfLen;i--,j++){
		wArr[i]=wArr[j];
	}

	return wArr;
}

static float *__window_createTukey(int length,float *alpha){
	float *wArr=NULL;
	float *xArr=NULL;

	float a=0.5;

	if(alpha){
		if(*alpha>=0&&*alpha<=1){
			a=*alpha;
		}
	}

	if(a==0){ // rect
		float _v=1;

		wArr=__vnew(length, &_v);
		return wArr;
	}
	else if(a==1){ // hann
		wArr=window_createHann(length, 0); // symmetric
		return wArr;
	}

	wArr=__vnew(length, NULL);
	xArr=__vlinspace(0, 1, length, 0);
	for(int i=0;i<length;i++){
		float _x=0;

		_x=xArr[i];
		if(_x>=0&&_x<a/2){
			wArr[i]=0.5*(1+cosf(2*M_PI/a*(_x-a/2)));
		}
		else if(_x>=a/2&&_x<(1-a/2)){
			wArr[i]=1;
		}
		else{ // 1-r/2~1
			wArr[i]=0.5*(1+cosf(2*M_PI/a*(_x-1+a/2)));
		}
	}

	free(xArr);
	return wArr;
}

static float *__calHalfBartlett(int halfLen,int length){
	float *arr=NULL;

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=0;i<halfLen;i++){
		arr[i]=2.0*i/length;
	}

	return arr;
}

// i=1 start
static float *__calHalfBartlettHann(int halfLen,int length){
	float *arr=NULL;

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=1;i<halfLen;i++){
		arr[i]=0.62-0.48*fabs(1.0*i/length-0.5)+0.38*cosf(2*M_PI*(1.0*i/length-0.5));
	}

	return arr;
}

static float *__calHalfTriang(int halfLen,int length){
	float *arr=NULL;

	int len=0;
	float det=0;

	if(!(length&1)){ // even
		det=0.5;
	}
	else{ // odd
		det=1;
		len=1;
	}

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=0;i<=halfLen;i++){
		arr[i]=2.0*(i+det)/(length+len);
	}

	return arr;
}

/****
	I(a)为零阶第一类修正贝塞尔函数
	I(a)=1+累加(k=1,无情大)[(1/k!(a/2)^k)^2]
	a一般取4~9,可以更小
****/
static float __besselZeroOne(float a){
	float sum=0;
	float b=0;

	float num=0;
	float den=0;
	float mid=0;

	sum=1;
	b=a/2;
	num=1;
	den=1;

	for(int k=1;k<16;k++){ // 一般取15~25 8即可满足精度
		num=num*b;
		den=den*k;
		mid=num/den;
		sum=sum+mid*mid;
	}

	return sum;
}

/***
	w(n)=I(b)/I(a) n>=0&&n<N-1;
	b=sqrt(1-(2*n/(N-1)-1)^2);
	a一般取4~9,可以更小
****/
static float *__calHalfKaiser(int halfLen,float a,int length){
	float *arr=NULL;

	float num=0;
	float den=0;
	float b=0;

	arr=__vnew(halfLen*2+1, NULL);
	den=__besselZeroOne(a);
	for(int i=0;i<halfLen;i++){
		float _value=0;

		_value=2.0*i/length-1;
		b=a*sqrtf(1-_value*_value);
		num=__besselZeroOne(b);
		arr[i]=num/den;
	}

	return arr;
}

// gauss函数
static float *__calHalfGauss(int halfLen,float a,int length){
	float *arr=NULL;
	float det=0;

	if(!(length&1)){ // even
		det=0.5;
	}

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=halfLen-1,j=0;i>=0;i--,j++){
		float _value=0;

		_value=2*a*(i-det)/(length-1);
		_value=-0.5*_value*_value;
		arr[j]=expf(_value);
	}

	return arr;
}

static float *__calHalfHann(int halfLen,int length){
	float *arr=NULL;

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=0;i<halfLen;i++){
		arr[i]=0.5-0.5*cosf(2*M_PI*i/length);
	}

	return arr;
}

static float *__calHalfHamm(int halfLen,int length){
	float *arr=NULL;

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=0;i<halfLen;i++){
		arr[i]=0.54-0.46*cosf(2*M_PI*i/length);
	}

	return arr;
}

// i=1 start
static float *__calHalfBlackman(int halfLen,int length){
	float *arr=NULL;

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=1;i<halfLen;i++){
		arr[i]=0.42-0.5*cosf(2*M_PI*i/length)+0.08*cosf(4*M_PI*i/length);
	}

	return arr;
}

// i=1 start
static float *__calHalfBohman(int halfLen,int length){
	float *arr=NULL;
	float *lArr=NULL;

	arr=__vnew(halfLen*2+1, NULL);
	lArr=__vlinspace(-1, 1, length, 0);
	for(int i=1;i<halfLen;i++){
		arr[i]=(1-fabsf(lArr[i]))*cosf(M_PI*fabsf(lArr[i]))+1/M_PI*sinf(M_PI*fabsf(lArr[i]));
	}

	free(lArr);
	return arr;
}

static float *__calHalfBlackmanHarris(int halfLen,int length){
	float *arr=NULL;

	float a0=0;
	float a1=0;
	float a2=0;
	float a3=0;

	a0=0.35875;
	a1=0.48829;
	a2=0.14128;
	a3=0.01168;

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=0;i<halfLen;i++){
		arr[i]=a0-a1*cosf(2*M_PI*i/length)+
				a2*cosf(4*M_PI*i/length)-
				a3*cosf(6*M_PI*i/length);
	}

	return arr;
}

static float *__calHalfBlackmanNuttall(int halfLen,int length){
	float *arr=NULL;

	float a0=0;
	float a1=0;
	float a2=0;
	float a3=0;

	a0=0.3635819;
	a1=0.4891775;
	a2=0.1365995;
	a3=0.0106411;

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=0;i<halfLen;i++){
		arr[i]=a0-a1*cosf(2*M_PI*i/length)+
				a2*cosf(4*M_PI*i/length)-
				a3*cosf(6*M_PI*i/length);
	}

	return arr;
}

static float *__calHalfFlattop(int halfLen,int length){
	float *arr=NULL;

	float a0=0;
	float a1=0;
	float a2=0;
	float a3=0;
	float a4=0;

	a0=0.21557895;
	a1=0.41663158;
	a2=0.277263158;
	a3=0.083578947;
	a4=0.006947368;

	arr=__vnew(halfLen*2+1, NULL);
	for(int i=0;i<halfLen;i++){
		arr[i]=a0-a1*cosf(2*M_PI*i/length)+
				a2*cosf(4*M_PI*i/length)-
				a3*cosf(6*M_PI*i/length)+
				a4*cosf(8*M_PI*i/length);
	}

	return arr;
}

void window_calKaiserOrder(float w1,float w2,float atten,int *order,float *beta){

	__calKaiserOrder(w1,w2,atten,order,beta);
}

static void __calKaiserOrder(float w1,float w2,float atten,int *order,float *beta){
	int _order=0;
	float _beta=0;

	if(w1<=0||w1>=1||w2<=0||w2>=1){
		return;
	}

	_order=ceilf((atten-7.95)/(M_PI*2.285*fabsf(w1-w2)));
	if(atten>50){
		_beta=0.1102*(atten-8.7);
	}
	else if(atten>=21){
		_beta=0.5842*powf((atten-21), 0.4)+
    			0.07886*(atten-21);
	}

	if(order){
		*order=_order;
	}

	if(beta){
		*beta=_beta;
	}
}

float *window_calFFTWindow(WindowType winType,int length){
	float *dataArr=NULL;

	if(winType==Window_Hann){
		dataArr=window_createHann(length,1);
	}
	else if(winType==Window_Hamm){
		dataArr=window_createHamm(length,1);
	}
	else if(winType==Window_Blackman){
		dataArr=window_createBlackman(length,1);
	}
	else if(winType==Window_Kaiser){ // beta 5
		dataArr=window_createKaiser(length,1,NULL);
	}
	else if(winType==Window_Bartlett){ // symmetric
		dataArr=window_createBartlett(length,0);
	}
	else if(winType==Window_Triang){ // symmetric
		dataArr=window_createTriang(length,0);
	}
	else if(winType==Window_Flattop){
		dataArr=window_createFlattop(length,1);
	}
	else if(winType==Window_Gauss){ // alpha 2.5
		dataArr=window_createGauss(length,1,NULL);
	}
	else if(winType==Window_Blackman_Harris){
		dataArr=window_createBlackmanHarris(length,1);
	}
	else if(winType==Window_Blackman_Nuttall){
		dataArr=window_createBlackmanNuttall(length,1);
	}
	else if(winType==Window_Bartlett_Hann){ // symmetric
		dataArr=window_createBartlettHann(length,0);
	}
	else if(winType==Window_Bohman){ // symmetric
		dataArr=window_createBohman(length,0);
	}
	else if(winType==Window_Tukey){ // alpha 0.5
		dataArr=window_createTukey(length,1,NULL);
	}
	else { // rect
		dataArr=(float *)calloc(length, sizeof(float ));
		for(int i=0;i<length;i++){
			dataArr[i]=1;
		}
	}

	return dataArr;
}






















