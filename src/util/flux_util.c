// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "../dsp/filterDesign_fir.h"

#include "../util/flux_wave.h"
#include "../util/flux_util.h"

static int __isEqual(float value1,float value2);
static int __isGreater(float value1,float value2);
static int __isLess(float value1,float value2);

// 2^n相关
int util_isPowerTwo(int value){
	int flag=1;

	if(value<1){
		flag=0;
	}

	if(value&(value-1)){
		flag=0;
	}

	return flag;
}

int util_ceilPowerTwo(int value){
	int n=1;

	if(value<1){
		return n;
	}

	if(util_isPowerTwo(value)){
		n=value;
	}
	else{
		while(value){
			value>>=1;
			n<<=1;
		}
	}

	return n;
}

int util_floorPowerTwo(int value){
	int n=1;

	if(value<1){
		return n;
	}

	if(util_isPowerTwo(value)){
		n=value;
	}
	else{
		value>>=1;
		while(value){
			value>>=1;
			n<<=1;
		}
	}

	return n;
}

int util_roundPowerTwo(int value){
	int n=1;

	int v1=1;
	int v2=1;

	if(value<1){
		return n;
	}
	
	if(util_isPowerTwo(value)){
		n=value;
	}
	else{
		v1=util_floorPowerTwo(value);
		v2=util_ceilPowerTwo(value);

		if(value-v1<v2-value){
			n=v1;
		}
		else{
			n=v2;
		}
	}

	return n;
}

// value=2^n =>n
int util_powerTwoBit(int value){
	int r=0;
	int m=0;

	if(util_isPowerTwo(value)){
		while(1){
			m=value%2;
			if(!m){ //
				r++;
				value=value/2;
				if(value==1){
					break;
				}
			}
			else{
				r=0;
				break;
			}
		}
	}

	return r;
}

// 最大公约数 a>b 辗转相除
int util_gcd(int a,int b){
	int c=0;

	c=a%b;
	if(c==0){
		return b;
	}
	else{
		return util_gcd(b,c);
	}
}

// tone
float util_freToMidiRange(float fre,int *midi1,int *midi2){
	float per=0;

	float curFre=0;
	float preFre=0;
	float nextFre=0;

	float sub1=0;
	float sub2=0;

	int m1=0;
	int m2=0;

	m1=round(12*log2(fre/440)+69);

	curFre=powf(2,1.0*(m1-69)/12)*440;
	preFre=powf(2,1.0*(m1-1-69)/12)*440;
	nextFre=powf(2,1.0*(m1+1-69)/12)*440;

	sub1=curFre-preFre;
	sub2=nextFre-curFre;

	if(nextFre-fre<fre-preFre){
		m2=m1+1;
	}
	else{
		m2=m1-1;
	}

	if(fre>curFre){
		if(m2==m1+1){
			per=(nextFre-fre)/sub2;
		}
		else{
			per=(sub1-(fre-curFre))/sub1;
		}
	}
	else{
		per=(sub1-(curFre-fre))/sub1;
	}

	if(midi1){
		*midi1=m1;
	}

	if(midi2){
		*midi2=m2;
	}

	return per;
}

void util_calTone(float value,float *value1,float *value2){
	float curFre=0;
	float preFre=0;
	float nextFre=0;

	int midi=0;

	midi=round(12*log2(value/440)+69);

	curFre=powf(2,1.0*(midi-69)/12)*440;
	if(value1){
		*value1=curFre;
	}

	if(value2){
		preFre=powf(2,1.0*(midi-1-69)/12)*440;
		nextFre=powf(2,1.0*(midi+1-69)/12)*440;

		if(nextFre-value<value-preFre){
			*value2=nextFre;
		}
		else{
			*value2=preFre;
		}
	}
}

int util_calToneTimes(float value1,float value2,int *type){
	int k=0;

	float _value=0;

	if(!value1||!value2){
		return 0;
	}

	if(type){
		*type=0;
	}

	if(__isEqual(value1, value2)){ // value1==value2
		k=1;
	}
	else if(__isLess(value1, value2)){ // value1<value2
		k=roundf(value2/value1);
		util_calTone(k*value1, &_value, NULL);
		if(!__isEqual(value2, _value)){
			k=0;
		}
	}
	else { // value1>value2
		k=roundf(value1/value2);
		util_calTone(k*value2, &_value, NULL);
		if(!__isEqual(value1, _value)){
			k=0;
		}

		if(type){
			*type=1;
		}
	}

	return k;
}

int util_calFreTimes(float value1,float value2,int *type){
	int flag=0;

	float _value1=0;
	float _value2=0;

	util_calTone(value1,&_value1,NULL);
	util_calTone(value2,&_value2,NULL);

	flag=util_calToneTimes(_value1,_value2,type);

	return flag;
}

/***
	1. sub:value
	2. sub1:sub2
****/
int util_calRangeTimes(float value1,float value2,int *type){
	int k=0;

	float _value1=0;
	float _value2=0;

	float _select1=0;
	float _select2=0;

	int flag1=0;
	int flag2=0;

	float s1=0;
	float s2=0;

	util_calTone(value1, &_value1, &_select1);
	util_calTone(value2, &_value2, &_select2);
	
	if(value1>660){
		s1=10;
	}
	else if(value1>330){
		s1=5;
	}

	if(fabsf(fabsf(_value1-value1)-fabsf(_select1-value1))<s1){
		flag1=1;
	}

	if(value2>660){
		s2=10;
	}
	else if(value2>330){
		s2=5;
	}

	if(fabsf(fabsf(_value2-value2)-fabsf(_select2-value2))<s2){
		flag2=1;
	}

	k=util_calToneTimes(_value1, _value2, type);
	if(!k&&(value1<330||flag1)){ // 338/340 330->440
		k=util_calToneTimes(_select1, _value2, type);
		if(!k&&(value2<330||flag2)){ // 338/340 330->440
			k=util_calToneTimes(_value1, _select2, type);
			if(!k){
				k=util_calToneTimes(_select1, _select2, type);
			}
		}
	}

	if(k>10){
		float e1=0,e2=0,e3=0;

		e1=fabsf((k-1)*value1-value2);
		e2=fabsf(k*value1-value2);
		e3=fabsf((k+1)*value1-value2);

		if(e1<e2&&e1<e3){
			k--;
		}
		else if(e3<e1&&e3<e2){
			k++;
		}
	}

	return k;
}

int util_calApproTimes(float value1,float value2,int *type){
	int k=0;

	return k;
}

// midi 12*log2(fre/440)+69
float util_midiToFre(int midi){
	float fre=0;

	fre=powf(2,1.0*(midi-69)/12)*440;

	return fre;
}

int util_freToMidi(float fre){
	int midi=0;

	midi=roundf(12*log2(fre/440)+69);

	return midi;
}

int util_midiTimes(int midi1,int midi2){
	int k=0;

	float fre1=0;
	float fre2=0;
	float fre3=0;

	int _m1=0;
	int _m2=0;

	if(midi1>=midi2){
		fre1=util_midiToFre(midi1);
		fre2=util_midiToFre(midi2);
		_m1=midi1;
	}
	else{
		fre1=util_midiToFre(midi2);
		fre2=util_midiToFre(midi1);
		_m1=midi2;
	}

	k=roundf(fre1/fre2);
	fre3=fre2*k;
	_m2=util_freToMidi(fre3);

	if(_m1!=_m2){
		k=0;
	}

	return k;
}

int util_freToSimularMidi(float fre){
	int midi1=0;
	int midi2=0;

	float tone1=0;
	float tone2=0;

	float det=0;
	float _m=0;

	midi1=util_freToMidi(fre);
	tone1=util_midiToFre(midi1);
	if(fre<tone1){
		midi2=midi1-1;
	}
	else{
		midi2=midi1+1;
	}

	tone2=util_midiToFre(midi2);
	det=tone1-tone2;
	_m=tone2+det/2;
	if(fabsf(fre-_m)>fabsf(det)/4){
		midi2=0;
	}

	return midi2;
}

int util_freTimes(float fre1,float fre2){
	int k=0;

	int midi1=0;
	int midi2=0;

	int _midi1=0;
	int _midi2=0;

	midi1=util_freToMidi(fre1);
	midi2=util_freToMidi(fre2);

	_midi1=util_freToSimularMidi(fre1);
	_midi2=util_freToSimularMidi(fre2);

	k=util_midiTimes(midi1,midi2);
	if(!k){
		if(midi1<midi2){
			if(_midi1){ // midi1<64
				k=util_midiTimes(_midi1,midi2);
			}

			if(!k&&_midi2){
				k=util_midiTimes(midi1,_midi2);
			}

			if(!k&&_midi1&&_midi2){
				k=util_midiTimes(_midi1,_midi2);
			}
		}
		else{
			if(_midi2){ // midi2<64
				k=util_midiTimes(midi1,_midi2);
			}

			if(!k&&_midi1){
				k=util_midiTimes(_midi1,midi2);
			}

			if(!k&&_midi1&&_midi2){
				k=util_midiTimes(_midi1,_midi2);
			}
		}
	}

	return k;
}

// scale
void util_minMaxScale(float *vArr1,int length,float *vArr2){

	__vminmaxscale(vArr1,length,vArr2);
}

void util_standScale(float *vArr1,int length,int type,float *vArr2){

	__vstandscale(vArr1,length,type,vArr2);
}

void util_maxAbsScale(float *vArr1,int length,float *vArr2){

	__vmaxabsscale(vArr1,length,vArr2);
}

void util_robustScale(float *vArr1,int length,float *vArr2){

	__vrobustscale(vArr1,length,vArr2);
}

void util_centerScale(float *vArr1,int length,float *vArr2){

	__vcenterscale(vArr1,length,vArr2);
}

void util_meanScale(float *vArr1,int length,float *vArr2){

	__vmeanscale(vArr1,length,vArr2);
}

void util_arctanScale(float *vArr1,int length,float *vArr2){

	__varctanscale(vArr1,length,vArr2);
}

// 向量归一化 type 0 p 1/2/3... 1 Inf 2 -Inf
void util_normalize(float *vArr1,int length,int type,float p,float *vArr2){

	__vnormalize(vArr1,length,type,p,vArr2);
}

float util_rectScale(float cur,float left,float right){
	float d1=0;
	float d2=0;

	d1=right-left;
	d2=(d1>=0)?(right)/(cur+right):(-left)/(cur+left);

	return d2;
}

float util_hannScale(float cur,float left,float right){
	float d1=0;
	float d2=0;

	d1=right-left;
	d2=(d1>=0)?(2*right-cur)/(cur+right):(cur-2*left)/(cur+left);

	return d2;
}

float util_hammScale(float cur,float left,float right){
	float d1=0;
	float d2=0;

	d1=right-left;
	d2=(d1>=0)?(2*right-cur)/(cur+right):(cur-2*left)/(cur+left);

	return d2;
}

void util_powerToDB(float *pArr,int length,float min,float *dArr){
	float max=0;
	float *arr=NULL;

	if(dArr){
		arr=dArr;
	}
	else{
		arr=pArr;
	}

	if(min>=0){
		min=-80;
	}

	__vmax(pArr, length, &max);
	for(int i=0;i<length;i++){
		float _value=0;

		_value=10*log10f(pArr[i]/max);
		arr[i]=(_value>min?_value:min);
	}
}

void util_powerToAbsDB(float *pArr,int length,int fftLength,int isNorm,float min,float *dArr){
	float max=0;
	int maxIndex=0;

	float *arr=NULL;
	float fLen=0;

	if(dArr){
		arr=dArr;
	}
	else{
		arr=pArr;
	}

	if(min>=0){
		min=-80;
	}

	fLen=fftLength*fftLength;
	for(int i=0;i<length;i++){
		float _value=0;

		_value=10*log10f(pArr[i]/fLen);
		arr[i]=(_value>min?_value:min);
	}

	if(isNorm){
		maxIndex=__vmax(pArr, length, NULL);
		max=arr[maxIndex];
		for(int i=0;i<length;i++){
			arr[i]=max-arr[i];
		}
	}
}

void util_magToAbsDB(float *pArr,int length,int fftLength,int isNorm,float min,float *dArr){
	float max=0;
	int maxIndex=0;

	float *arr=NULL;

	if(dArr){
		arr=dArr;
	}
	else{
		arr=pArr;
	}

	if(min>=0){
		min=-80;
	}

	for(int i=0;i<length;i++){
		float _value=0;

		_value=20*log10f(pArr[i]/fftLength);
		arr[i]=(_value>min?_value:min);
	}

	if(isNorm){
		maxIndex=__vmax(pArr, length, NULL);
		max=arr[maxIndex];
		for(int i=0;i<length;i++){
			arr[i]=max-arr[i];
		}
	}
}

void util_logCompress(float *vArr1,float *gamma,int length,float *vArr2){

	__vlog_compress(vArr1, gamma,NULL, length, vArr2);
}

void util_log10Compress(float *vArr1,float *gamma,int length,float *vArr2){

	__vlog10_compress(vArr1, gamma,NULL, length, vArr2);
}

// temproal, return maxDb, percent is < -base
float util_temproal(float *dataArr,int length,float base,float *avgValue,float *percentValue){
	float max=-100;

	int count=0;
	float value=0;
	float total=0;
	
	if(!dataArr||!length){
		return 0;
	}

	for(int i=0;i<length;i++){
		value=20*log10f(fabsf(dataArr[i])+1e-8);
		if(value<-36){
			value=-36;
		}

		if(max<value){
			max=value;
		}

		total+=value;

		if(value>-base){
			count++;
		}
	}

	*avgValue=total/length;
	*percentValue=1.0*(length-count)/length;

	return max;
}

/***
	x=1 1
	0<x<1 
	x>1 递归
	x<0 递归
****/
float util_gamma(float x,float *eps){
	float y=0;
	float err=0;

	int i=1;
	float cur=1.0;
	float pre=0;

	err=1e-5;
	if(eps){
		if(*eps>0){
			err=*eps;
		}
	}

	if(fabs(x-1.0)<err){
		return 1.0;
	}
	else if(fabs(x-0.5)<err){
		return sqrt(M_PI);
	}

	if(x>1.0){
		return (x-1)*util_gamma(x-1, eps);
	}
	else if(x<0){
		return util_gamma(x+1, eps)/x;
	}

	while(fabsf(cur-pre)/cur>err){
		pre=cur;
		cur*=i/(x-1+i);
		i++;
	}
	y=cur*powf(i, x-1);

	return y;
}

long double util_gammal(long double x,long double *eps){
	long double y=0;
	long double err=0;

	int i=1;
	long double cur=1.0;
	long double pre=0;

	err=1e-5;
	if(eps){
		if(*eps>0){
			err=*eps;
		}
	}

	if(fabsl(x-1.0)<err){
		return 1.0;
	}
	else if(fabsl(x-0.5)<err){
		return sqrt(M_PI);
	}

	if(x>1.0){
		return (x-1)*util_gammal(x-1, eps);
	}
	else if(x<0){
		return util_gammal(x+1, eps)/x;
	}

	while(fabsl(cur-pre)/cur>err){
		pre=cur;
		cur*=i/(x-1+i);
		i++;
	}
	y=cur*powl(i, x-1);

	return y;
}

// quadratically interpolated return p ∈[-1/2,1/2]
float util_qaudInterp(float value1,float value2,float value3,float *outValue){
	float p=0;

	p=(value3-value1)/(2*(2*value2-value3-value1)+1e-16);
	if(outValue){
		*outValue=value2-0.25*(value1-value3)*p;
	}

	return p;
}

// peak pick
void util_peakPick(float *pArr,int length,int start,int end,int distance,int num,float *vArr,int *cArr){
	int index1=0;
	int index2=0;

	if(distance<=0){
		distance=1;
	}

	for(int i=0;i<num;i++){
		cArr[i]=__vmax(pArr+start, end-start+1, vArr+i)+start;

		index1=(cArr[i]-distance>start?cArr[i]-distance:start);
		index2=(cArr[i]+distance<end?cArr[i]+distance:end);
		for(int j=index1;j<=index2;j++){
			pArr[j]=NAN;
		}
	}
}

// order must odd; delta/deltaDelta
void util_delta(float *dataArr1,int length,int order,float *dataArr2){
	float *bArr=NULL;
	float aArr[]={1};

	if(order&1){
		bArr=filterDesign_smooth1(order);
		filterDesign_filter(bArr,aArr,dataArr1,
							order,1,length,
							dataArr2);
	}

	free(bArr);
}

// pre emphasis coef 0.97
void util_preEmphasis(float *vArr1,int length,float coef,float *vArr2){

	if(vArr2){
		vArr2[0]=vArr1[0];
		for(int i=1;i<length;i++){
			vArr2[i]=vArr1[i]-coef*vArr1[i-1];
		}
	}
}

int util_readWave(char *name,float **dataArr){
	int status=0;
	WaveReadObj waveRead=NULL;

	int samplate=0;
	int bit=0;
	int channelNum=0;

	int length=0;
	float *_dataArr=NULL;
	
	status=waveReadObj_new(&waveRead,name);
	if(!status){
		length=waveReadObj_getInfor(waveRead,&samplate,&bit,&channelNum);
		_dataArr=__vnew(length, NULL);

		waveReadObj_read(waveRead, _dataArr, length);
		*dataArr=_dataArr;

		waveReadObj_free(waveRead);
	}

	return length;
}

void util_writeWave(char *name,float *dataArr,int length){
	int status=0;

	WaveWriteObj waveWrite=NULL;

	int samplate=32000;
	int bit=16;
	int channelNum=1;

	status=waveWriteObj_new(&waveWrite,name,
							&samplate,&bit,&channelNum);
	if(!status){
		waveWriteObj_write(waveWrite,dataArr,length);
		waveWriteObj_free(waveWrite);
	}
}

static int __isEqual(float value1,float value2){
	int flag=0;

	if(fabsf(value1-value2)<0.81){
		return 1;
	}

	return flag;
}

static int __isGreater(float value1,float value2){
	int flag=0;

	if(value1-value2>0.81){
		return 1;
	}

	return flag;
}

static int __isLess(float value1,float value2){
	int flag=0;

	if(value2-value1>0.81){
		return 1;
	}


	return flag;
}









