// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../dsp/flux_window.h"
#include "../dsp/filterDesign_freqz.h"

#include "auditory_filterBank.h"

// ESTI
static void __auditory_estiFilterBank(int num,int fftLength,int isPseudo,
									SpectralFilterBankNormalType type,
									float *freBandArr,int *binBandArr,
									float *mFilterBankArr);
// Slaney
static void __auditory_slaneyFilterBank(int num,int fftLength,int samplate,int isPseudo,
									SpectralFilterBankNormalType type,
									float *freBandArr,int *binBandArr,
									float *mFilterBankArr);
// gammatone
static void __auditory_gammatoneFilterBank(int num,int fftLength,int samplate,int isPseudo,
									SpectralFilterBankNormalType type,
									float *freBandArr,int *binBandArr,
									float *mFilterBankArr);

// rect~guass 窗函数设计法
static void __auditory_windowFilterBank(int num,int fftLength,int samplate,int isPseudo,
									SpectralFilterBankStyleType styleType,
									SpectralFilterBankNormalType normType,
									float *freBandArr,int *binBandArr,
									float *mFilterBankArr);

// linear
static void __auditory_linearFilterBank(int num,int fftLength,int samplate,int isPseudo,
										float *freBandArr,int *binBandArr,
										float *mFilterBankArr);

// isEdge 0 针对非gammatone 1 gammatone
static void __auditory_calBandEdge(int num,int fftLength,int samplate,
								float lowFre,float highFre,int isEdge,int bType,
								void *callback1,void *callback2,float ref,
								float **freBandArr,int **binBandArr);

static void __reviseMidiFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3);
static void __reviseLogFre(int num,float lowFre,float highFre,int binPerOctave,int isEdge,float *lowFre3,float *highFre3);
static void __reviseLinearFre(int num,float lowFre,float highFre,float detFre,int isEdge,float *lowFre3,float *highFre3);

static void __reviseLinspaceFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3);
static void __reviseLogspaceFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3);

void auditory_filterBank(int num,int fftLength,int samplate,int isPseudo,
						SpectralFilterBankScaleType scaleType,SpectralFilterBankStyleType styleType,SpectralFilterBankNormalType normType,
						float lowFre,float highFre,int binPerOctave,
						float *mFilterBankArr,
						float *freBandArr,
						int *binBandArr){
	float *fArr=NULL;
	int *bArr=NULL;

	int isEdge=0;
	int offset=0;

	// UniFunc||UniFunc1
	void *func1=NULL;
	void *func2=NULL;

	float ref=0;
	int bType=0;

	if(styleType==SpectralFilterBankStyle_Gammatone){ // include edge!!!
		isEdge=1;
	}

	if(!isEdge){
		offset=1;
	}

	// 0. revise log/linear!!!
	if(scaleType==SpectralFilterBankScale_Octave){ 
		if(binPerOctave&&binPerOctave>=4&&binPerOctave<=48){
			ref=binPerOctave;
		}
		else{
			ref=12;
		}

		__reviseLogFre(num,lowFre,highFre,ref,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Linear){
		ref=samplate*1.0/fftLength;

		__reviseLinearFre(num,lowFre,highFre,ref,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Linspace){
		__reviseLinspaceFre(num,lowFre,highFre,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Log){
		__reviseLogspaceFre(num,lowFre,highFre,isEdge,&lowFre,&highFre);
	}

	// compatible LogChroma !!!
	if(scaleType==SpectralFilterBankScale_LogChroma){ 
		if(binPerOctave>=12&&binPerOctave%12==0){
			ref=binPerOctave;
		}
		else{
			ref=12;
		}

		__reviseLogFre(num,lowFre,highFre,ref,isEdge,&lowFre,&highFre);
	}

	// 1. freBand/binBand
	if(scaleType==SpectralFilterBankScale_Linear){ // linear
		func1=auditory_freToLinear;
		func2=auditory_linearToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Linspace){ // linspace
		func1=auditory_freToLinspace;
		func2=auditory_linspaceToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Mel){ // mel
		func1=auditory_freToMel;
		func2=auditory_melToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Bark){ // bark
		func1=auditory_freToBark;
		func2=auditory_barkToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Erb){ // erb
		func1=auditory_freToErb;
		func2=auditory_erbToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Octave){ // log
		func1=auditory_freToLog;
		func2=auditory_logToFre;
	}
	else if(scaleType==SpectralFilterBankScale_Log){ // logspace
		func1=auditory_freToLogspace;
		func2=auditory_logspaceToFre;
	}
	else if(scaleType==SpectralFilterBankScale_LogChroma){ // compatible LogChroma !!!
		func1=auditory_freToLog;
		func2=auditory_logToFre;
	}

	if(styleType==SpectralFilterBankStyle_Slaney){
		bType=1;
	}

	__auditory_calBandEdge(num,fftLength,samplate,
						lowFre,highFre,isEdge,bType,
						func1,func2,ref,
						&fArr,&bArr);
	
	// 2. filterBank
	if(scaleType==SpectralFilterBankScale_Linear){
		__auditory_linearFilterBank(num,fftLength,samplate,isPseudo,
									fArr,bArr,
									mFilterBankArr);
	}
	else{
		if(styleType==SpectralFilterBankStyle_Slaney){ // Slaney style

			__auditory_slaneyFilterBank(num,fftLength,samplate,isPseudo,
										normType,
										fArr,bArr,
										mFilterBankArr);
		}
		else if(styleType==SpectralFilterBankStyle_ETSI){ // ESTI style

			__auditory_estiFilterBank(num,fftLength,isPseudo,
									normType,
									fArr,bArr,
									mFilterBankArr);
		}
		else if(styleType==SpectralFilterBankStyle_Gammatone){ // gammatone style
			__auditory_gammatoneFilterBank(num,fftLength,samplate,isPseudo,
										normType,
										fArr,bArr,
										mFilterBankArr);
		}
		else{ // rect~gauss
			__auditory_windowFilterBank(num,fftLength,samplate,isPseudo,
										styleType,normType,
										fArr,bArr,
										mFilterBankArr);
		}
	}

	// 3. 处理freBand/binBand
	if(freBandArr){
		memcpy(freBandArr, fArr+offset, sizeof(float )*num);
	}
	
	if(binBandArr){
		memcpy(binBandArr, bArr+offset, sizeof(int )*num);
	}

	free(fArr);
	free(bArr);
}

// rect~guass 窗函数设计法
static void __auditory_windowFilterBank(int num,int fftLength,int samplate,int isPseudo,
									SpectralFilterBankStyleType winType,
									SpectralFilterBankNormalType normType,
									float *freBandArr,int *binBandArr,
									float *mFilterBankArr){
	int left=0;
	int cur=0;
	int right=0;

	int mLength=0;

	if(!isPseudo){
		mLength=fftLength/2+1;
	}
	else{
		mLength=fftLength;
	}

	// 1. all window freband
	if(winType==SpectralFilterBankStyle_Point){
		for(int i=1;i<num+1;i++){
			left=binBandArr[i-1];
			cur=binBandArr[i];
			right=binBandArr[i+1];

			mFilterBankArr[(i-1)*mLength+cur]=1.0;
		}
	}
	else if(winType==SpectralFilterBankStyle_Rect){
		for(int i=1;i<num+1;i++){
			left=binBandArr[i-1];
			cur=binBandArr[i];
			right=binBandArr[i+1];

			for(int j=left;j<=right;j++){
				mFilterBankArr[(i-1)*mLength+j]=1.0;
			}
		}
	}
	else{
		for(int i=1;i<num+1;i++){
			left=binBandArr[i-1];
			cur=binBandArr[i];
			right=binBandArr[i+1];

			if(cur>left){
				float *wArr=NULL;
				int _index=0;

				if(winType==SpectralFilterBankStyle_Hann){
					wArr=window_createHann(2*(cur-left)+1,0);
				}
				else if(winType==SpectralFilterBankStyle_Hamm){
					wArr=window_createHamm(2*(cur-left)+1,0);
				}
				else if(winType==SpectralFilterBankStyle_Blackman){
					wArr=window_createBlackman(2*(cur-left)+1,0);
				}
				else if(winType==SpectralFilterBankStyle_Bohman){
					wArr=window_createBohman(2*(cur-left)+1,0);
				}
				else if(winType==SpectralFilterBankStyle_Kaiser){
					wArr=window_createKaiser(2*(cur-left)+1,0,NULL);
				}
				else{ // gauss
					wArr=window_createGauss(2*(cur-left)+1,0,NULL);
				}

				for(int j=left,k=_index;j<=cur;j++,k++){ 
					mFilterBankArr[(i-1)*mLength+j]=wArr[k];
				}

				free(wArr);
			}

			if(right>cur){
				float *wArr=NULL;
				int _index=0;

				_index=(2*(right-cur)+1)/2+1;
				if(winType==SpectralFilterBankStyle_Hann){
					wArr=window_createHann(2*(right-cur)+1,0);
				}
				else if(winType==SpectralFilterBankStyle_Hamm){
					wArr=window_createHamm(2*(right-cur)+1,0);
				}
				else if(winType==SpectralFilterBankStyle_Blackman){
					wArr=window_createBlackman(2*(right-cur)+1,0);
				}
				else if(winType==SpectralFilterBankStyle_Bohman){
					wArr=window_createBohman(2*(right-cur)+1,0);
				}
				else if(winType==SpectralFilterBankStyle_Kaiser){
					wArr=window_createKaiser(2*(right-cur)+1,0,NULL);
				}
				else{ // gauss
					wArr=window_createGauss(2*(right-cur)+1,0,NULL);
				}
				
				for(int j=cur+1,k=_index;j<=right;j++,k++){ 
					mFilterBankArr[(i-1)*mLength+j]=wArr[k];
				}

				free(wArr);
			}
		}
	}

	// 2. 归一化
	if(normType==SpectralFilterBankNormal_Area||normType==SpectralFilterBankNormal_BandWidth){
		float *weightArr=NULL;

		weightArr=__vnew(num, NULL);
		if(normType==SpectralFilterBankNormal_Area){ // area
			__msum(mFilterBankArr, num, mLength, 1, weightArr);
			// mel/bark 计算一半
			// __vdiv_value(weightArr, 2, num, NULL);
		}
		else{ // bandwidth
			__vsub(freBandArr+2, freBandArr, num, weightArr);
			__vdiv_value(weightArr, 2, num, NULL);
		}

		__mdiv_vector(mFilterBankArr, weightArr,0, num, mLength, 1, NULL);

		free(weightArr);
	}
}

static void __auditory_linearFilterBank(int num,int fftLength,int samplate,int isPseudo,
										float *freBandArr,int *binBandArr,
										float *mFilterBankArr){
	int cur=0;
	int mLength=0;

	if(!isPseudo){
		mLength=fftLength/2+1;
	}
	else{
		mLength=fftLength;
	}

	// ???
	// for(int i=0;i<num;i++){
	// 	cur=binBandArr[i];
	// 	mFilterBankArr[i*mLength+cur]=1;
	// }

	// linear spec,binBandArr-1
	for(int i=1;i<num+1;i++){
		binBandArr[i]-=1;

		cur=binBandArr[i];
		mFilterBankArr[(i-1)*mLength+cur]=1;
	}
}

/***
	针对mel/bark ESTI-style 'bin'
	1. 三角频带
	2. 归一化
	filterBankArr=>matrix num*(fftLength/2+1)
****/
static void __auditory_estiFilterBank(int num,int fftLength,int isPseudo,
									SpectralFilterBankNormalType type,
									float *freBandArr,int *binBandArr,
									float *mFilterBankArr){
	int left=0;
	int cur=0;
	int right=0;

	int mLength=0;

	if(!isPseudo){
		mLength=fftLength/2+1;
	}
	else{
		mLength=fftLength;
	}

	// 1. 三角频带
	for(int i=1;i<num+1;i++){
		left=binBandArr[i-1];
		cur=binBandArr[i];
		right=binBandArr[i+1];

		if(cur>left){
			for(int j=left;j<=cur;j++){ // 三角上升 
				mFilterBankArr[(i-1)*mLength+j]=1.0*(j-left)/(cur-left);
			}
		}

		for(int j=cur+1;j<=right;j++){ // 三角下降
			mFilterBankArr[(i-1)*mLength+j]=1.0*(right-j)/(right-cur);
		}
	}

	// 2. norm
	if(type==SpectralFilterBankNormal_Area||type==SpectralFilterBankNormal_BandWidth){
		float *weightArr=NULL;

		weightArr=__vnew(num, NULL);
		if(type==SpectralFilterBankNormal_Area){ // area
			__msum(mFilterBankArr, num, mLength, 1, weightArr);
			// mel/bark 计算一半
			// __vdiv_value(weightArr, 2, num, NULL);
		}
		else{ // bandwidth
			__vsub(freBandArr+2, freBandArr, num, weightArr);
			__vdiv_value(weightArr, 2, num, NULL);
		}

		__mdiv_vector(mFilterBankArr, weightArr,0, num, mLength, 1, NULL);

		free(weightArr);
	}
}

/***
	针对mel/bark Slaney-style 'fre'
	1. 确定顶点
	2. 三角频带
	3. 归一化
	filterBankArr=>matrix num*(fftLength/2+1)
****/
static void __auditory_slaneyFilterBank(int num,int fftLength,int samplate,int isPseudo,
									SpectralFilterBankNormalType type,
									float *freBandArr,int *binBandArr,
									float *mFilterBankArr){
	float *fArr=NULL;
	float *wArr=NULL;

	int mLength=0;

	if(!isPseudo){
		mLength=fftLength/2+1;
	}
	else{
		mLength=fftLength;
	}

	wArr=__vnew(num+1, NULL);

	// 1. linear freArr; not samplate/2.0,fftLength/2+1 for stop overflow
	fArr=__vlinspace(0, samplate-samplate/(float )fftLength, fftLength, 0);
	// fArr=__vlinspace(0, samplate/2.0, fftLength/2+1, 0);

	// for(int i=0;i<num+2;i++){
	// 	for(int j=0;j<fftLength;j++){
	// 		if(fArr[j]>freBandArr[i]){
	// 			bArr[i]=j;
	// 			break;
	// 		}
	// 	}
	// }

	// 2. slaney
	__vsub(freBandArr+1, freBandArr, num+1, wArr);
	for(int i=0;i<num;i++){
		for(int j=binBandArr[i];j<=binBandArr[i+1]-1;j++){ 
			mFilterBankArr[i*mLength+j]=(fArr[j]-freBandArr[i])/wArr[i];
		}

		for(int j=binBandArr[i+1];j<=binBandArr[i+2]-1;j++){ 
			mFilterBankArr[i*mLength+j]=(freBandArr[i+2]-fArr[j])/wArr[i+1];
		}
	}

	// 3. norm
	if(type==SpectralFilterBankNormal_Area||type==SpectralFilterBankNormal_BandWidth){
		float *weightArr=NULL;

		weightArr=__vnew(num, NULL);
		if(type==SpectralFilterBankNormal_Area){ // area
			__msum(mFilterBankArr, num, mLength, 1, weightArr);
			// mel/bark 计算一半
			// __vdiv_value(weightArr, 2, num, NULL);
		}
		else{ // bandwidth
			__vsub(freBandArr+2, freBandArr, num, weightArr);
			__vdiv_value(weightArr, 2, num, NULL);
		}

		__mdiv_vector(mFilterBankArr, weightArr,0, num, mLength, 1, NULL);

		free(weightArr);
	}

	free(fArr);
	free(wArr);
}

/***
	针对erb
	1. freBand->coeff 频带计算滤波器系数
	2. coeff->bank 
	3. 归一化
	4. scale down/up 0&&fftLength/2两端不变 中间*2
****/
static void __auditory_gammatoneFilterBank(int num,int fftLength,int samplate,int isPseudo,
										SpectralFilterBankNormalType type,
										float *freBandArr,int *binBandArr,
										float *mFilterBankArr){
	float **cArrArr=NULL;

	int mLength=0;

	float *rArr=NULL;
	float *iArr=NULL;

	int isWhole=0;
	int len=0;

	int order=4;

	if(!isPseudo){
		mLength=fftLength/2+1;
	}
	else{
		mLength=fftLength;
	}
	
	if(!isWhole){
		len=fftLength/2+1;
	}
	else{
		len=fftLength;
	}

	// 1. freBand ->coeff 
	cArrArr=auditory_calGammatoneCoefficient(freBandArr,num,samplate);

	// 2. coeff->bank
	__vcnew(len, NULL, &rArr, &iArr);
	for(int i=0;i<num;i++){ // 32
		filterDesign_freqzSOS(cArrArr[i], order, 
							fftLength, samplate, isWhole,
							NULL,
							rArr, iArr, NULL);

		__vcabs(rArr, iArr, len, mFilterBankArr+i*len);
	}

	// 3. 归一化
	if(type==SpectralFilterBankNormal_Area||type==SpectralFilterBankNormal_BandWidth){
		float *weightArr=NULL;
		
		weightArr=__vnew(num, NULL);
		if(type==SpectralFilterBankNormal_Area){ // area
			if(!isWhole){ // 0~fftLength/2+1
				for(int i=0;i<num;i++){ // 0+fftLength/2
					weightArr[i]=mFilterBankArr[i*mLength]+mFilterBankArr[i*mLength+mLength-1];
					weightArr[i]+=__vsum(mFilterBankArr+(i*mLength+1), mLength-2)*2;
				}
			}
			else{
				__msum(mFilterBankArr, num, mLength, 1, weightArr);
			}
		}
		else{ // bandwidth
			for(int i=0;i<num;i++){
				weightArr[i]=1.019*24.7*(0.00437*freBandArr[i]+1);
			}

			__vdiv_value(weightArr, 2, num, NULL);
		}

		__mdiv_vector(mFilterBankArr, weightArr,0, num, mLength, 1, NULL);

		free(weightArr);
	}

	// 4. scale down/up
	for(int i=0;i<num;i++){
		for(int j=1;j<mLength-1;j++){
			mFilterBankArr[i*mLength+j]*=2;
		}
	}

	free(rArr);
	free(iArr);
}

// 计算band edge 针对mel/bark/erb/log bType 0 !slaney 1 slaney
static void __auditory_calBandEdge(int num,int fftLength,int samplate,
								float lowFre,float highFre,int isEdge,int bType,
								void *callback1,void *callback2,float ref,
								float **freBandArr,int **binBandArr){
	float *fArr=NULL;
	int *bArr=NULL;

	float low=0;
	float high=0;

	int det=0;

	UniFunc func1=NULL;
	UniFunc func2=NULL;

	UniFunc1 func3=NULL;
	UniFunc1 func4=NULL;

	if(!isEdge){ // ! gammatone
		det=2;
	}

	if(ref){
		func3=(UniFunc1 )callback1;
		func4=(UniFunc1 )callback2;

		low=func3(lowFre,ref);
		high=func3(highFre,ref);
	}
	else{
		func1=(UniFunc )callback1;
		func2=(UniFunc )callback2;

		low=func1(lowFre);
		high=func1(highFre);
	}

	// 1. mel/barkArr
	fArr=__vlinspace(low, high, num+det, 0);

	// 2. mel/barkArr=>freArr
	if(ref){
		__vmap1(fArr, num+det, func4,ref, NULL);
	}
	else{
		__vmap(fArr, num+det, func2, NULL);
	}

	if(freBandArr){
		*freBandArr=fArr;
	}

	// 3. freArr=>binArr vector???
	if(binBandArr){
		bArr=__vnewi(num+det, NULL);
		if(!bType){ // 非slaney
			for(int i=0;i<num+det;i++){
				bArr[i]=roundf(fftLength*fArr[i]/samplate);
			}
		}
		else{ // slaney
			float *_arr1=NULL;

			// not samplate/2.0,fftLength/2+1 for stop overflow
			_arr1=__vlinspace(0, samplate-samplate/(float )fftLength, fftLength, 0);
			for(int i=0;i<num+2;i++){
				for(int j=0;j<fftLength;j++){
					if(_arr1[j]>fArr[i]){
						bArr[i]=j;
						break;
					}
				}
			}

			free(_arr1);
		}

		*binBandArr=bArr;
	}

	if(!freBandArr){
		free(fArr);
	}
}

/***
	gammatone filter
		g(t)=(a*t^(n-1))*(e^(-2*PI*b*t))*cos(2*PI*fc*t+p)
		a= amp factor
		n= filter order,一般选用4阶(人类听觉系统模型)
		fc= center fre
		b= bandwidth ,1.019*fre2Erb(fc)
		p= phase factor
	1. 计算A11 A23 A13 A14 arr;A0 A2 B0 B1 B2
	2. 计算gain Arr共10个num向量
	3. 有10个num Arr 拼接num个4*6 matrix
****/
float **auditory_calGammatoneCoefficient(float *freBandArr,int length,int samplate){
	float **coefArrArr=NULL;

	float *erbArr=NULL; // 1.019*fre2Erb(fre)
	float *argArr=NULL; // 2*PI*fre*t
	float *vArr=NULL; // -t*e^(-(erbArr*t))

	float *cosArr=NULL; // cos(argArr)
	float *sinArr=NULL; // sin(argArr)

	// 复数相关
	float *cRealArr=NULL; // -2*t*e^(i*2*argArr)
	float *cImageArr=NULL;
	float *gRealArr=NULL; // 2*t*e^(i*argArr-erbArr*t)
	float *gImageArr=NULL;

	float *k11Arr=NULL;
	float *k12Arr=NULL;
	float *k13Arr=NULL;
	float *k14Arr=NULL;

	// BA gain系数相关
	float *a0Arr=NULL; // 1/samplate
	float *a2Arr=NULL; // 0
	float *b0Arr=NULL; // 1

	float *b1Arr=NULL;
	float *b2Arr=NULL;

	float *a11Arr=NULL;
	float *a12Arr=NULL;
	float *a13Arr=NULL;
	float *a14Arr=NULL;

	float *gainArr=NULL;

	float t=0;
	float v1=1;

	float pv=0;
	float nv=0;

	int nLen=0;
	int mLen=0;

	t=1.0/samplate;
	pv=sqrtf(3+powf(2, 1.5));
	nv=sqrtf(3-powf(2, 1.5));

	nLen=4;
	mLen=6;

	coefArrArr=(float **)calloc(length, sizeof(float *));
	for(int i=0;i<length;i++){
		coefArrArr[i]=__vnew(nLen*mLen, NULL);
	}

	erbArr=__vnew(length, NULL);
	argArr=__vnew(length, NULL);
	vArr=__vnew(length, NULL);

	cosArr=__vnew(length, NULL);
	sinArr=__vnew(length, NULL);

	cRealArr=__vnew(length, NULL);
	cImageArr=__vnew(length, NULL);
	gRealArr=__vnew(length, NULL);
	gImageArr=__vnew(length, NULL);

	k11Arr=__vnew(length, &t);
	k12Arr=__vnew(length, &t);
	k13Arr=__vnew(length, &t);
	k14Arr=__vnew(length, &t);

	a0Arr=__vnew(length, &t);
	a2Arr=__vnew(length, NULL);
	b0Arr=__vnew(length, &v1);

	b1Arr=__vnew(length, NULL);
	b2Arr=__vnew(length, NULL);

	a11Arr=__vnew(length, NULL);
	a12Arr=__vnew(length, NULL);
	a13Arr=__vnew(length, NULL);
	a14Arr=__vnew(length, NULL);

	gainArr=__vnew(length, NULL);

	// erbArr argArr
	for(int i=0;i<length;i++){
		erbArr[i]=(freBandArr[i]/9.26449+24.7)*2*M_PI*1.019; // fre=>erb space
		argArr[i]=freBandArr[i]*2*M_PI*t;
	}

	// vArr
	for(int i=0;i<length;i++){
		vArr[i]=-t*expf(-t*erbArr[i]);
	}

	// cosArr/sinArr
	for(int i=0;i<length;i++){
		cosArr[i]=cosf(argArr[i]);
		sinArr[i]=sinf(argArr[i]);
	}

	// real/image arr
	for(int i=0;i<length;i++){
		cRealArr[i]=cosf(4*M_PI*t*freBandArr[i]);
		cImageArr[i]=sinf(4*M_PI*t*freBandArr[i]);

		gRealArr[i]=2*t*expf(-erbArr[i]*t)*cos(2*M_PI*t*freBandArr[i]);
		gImageArr[i]=2*t*expf(-erbArr[i]*t)*sin(2*M_PI*t*freBandArr[i]);
	}

	// b1/b2 arr
	for(int i=0;i<length;i++){
		b1Arr[i]=-2*cosArr[i]/expf(erbArr[i]*t);
		b2Arr[i]=expf(-2*t*erbArr[i]);
	}

	// a11/a12/a13/a14
	for(int i=0;i<length;i++){
		k11Arr[i]=cosArr[i]+pv*sinArr[i];
		k12Arr[i]=cosArr[i]-pv*sinArr[i];
		k13Arr[i]=cosArr[i]+nv*sinArr[i];
		k14Arr[i]=cosArr[i]-nv*sinArr[i];

		a11Arr[i]=vArr[i]*k11Arr[i];
		a12Arr[i]=vArr[i]*k12Arr[i];
		a13Arr[i]=vArr[i]*k13Arr[i];
		a14Arr[i]=vArr[i]*k14Arr[i];
	}

	// gainArr
	for(int i=0;i<length;i++){
		float value=0;

		float r1=0,i1=0;
		float r2=0,i2=0;
		float r3=0,i3=0;
		float r4=0,i4=0;
		float r5=0,i5=0;

		r1=-2*t*cRealArr[i]+gRealArr[i]*k11Arr[i];
		i1=-2*t*cImageArr[i]+gImageArr[i]*k11Arr[i];

		r2=-2*t*cRealArr[i]+gRealArr[i]*k12Arr[i];
		i2=-2*t*cImageArr[i]+gImageArr[i]*k12Arr[i];

		r3=-2*t*cRealArr[i]+gRealArr[i]*k13Arr[i];
		i3=-2*t*cImageArr[i]+gImageArr[i]*k13Arr[i];

		r4=-2*t*cRealArr[i]+gRealArr[i]*k14Arr[i];
		i4=-2*t*cImageArr[i]+gImageArr[i]*k14Arr[i];

		r5=-2/expf(2*t*erbArr[i])-2*cRealArr[i]+2*(1+cRealArr[i])/expf(t*erbArr[i]);
		i5=-2*cImageArr[i]+2*cImageArr[i]/expf(t*erbArr[i]);

		value=sqrtf(r1*r1+i1*i1)*
				sqrtf(r2*r2+i2*i2)*
				sqrtf(r3*r3+i3*i3)*
				sqrtf(r4*r4+i4*i4)/((r5*r5+i5*i5)*(r5*r5+i5*i5));
		// value/=__complexPowM(r5,i5,4);

		gainArr[i]=value;
	}

	// 拼接length个n*m matraix
	for(int k=0;k<length;k++){
		for(int i=0;i<nLen;i++){ // nLength
			if(i==0){
				coefArrArr[k][0]=a0Arr[k]/gainArr[k];
				coefArrArr[k][1]=a11Arr[k]/gainArr[k];
				coefArrArr[k][2]=a2Arr[k]/gainArr[k];
				coefArrArr[k][3]=b0Arr[k];
				coefArrArr[k][4]=b1Arr[k];
				coefArrArr[k][5]=b2Arr[k];
			}
			else{
				float *arr=NULL;

				if(i==1){
					arr=a12Arr;
				}
				else if(i==2){
					arr=a13Arr;
				}
				else {
					arr=a14Arr;
				}

				coefArrArr[k][i*mLen+0]=a0Arr[k];
				coefArrArr[k][i*mLen+1]=arr[k];
				coefArrArr[k][i*mLen+2]=a2Arr[k];
				coefArrArr[k][i*mLen+3]=b0Arr[k];
				coefArrArr[k][i*mLen+4]=b1Arr[k];
				coefArrArr[k][i*mLen+5]=b2Arr[k];
			}
		}
	}

	free(erbArr);
	free(argArr);
	free(vArr);

	free(cosArr);
	free(sinArr);

	free(cRealArr);
	free(cImageArr);
	free(gRealArr);
	free(gImageArr);

	free(k11Arr);
	free(k12Arr);
	free(k13Arr);
	free(k14Arr);

	free(a0Arr);
	free(a2Arr);
	free(b0Arr);

	free(b1Arr);
	free(b2Arr);

	free(a11Arr);
	free(a12Arr);
	free(a13Arr);
	free(a14Arr);

	free(gainArr);

	return coefArrArr;
}

// isEdge 0 非gammatone high=low+num-1+2; 1 gammatone high=low+num-1
static void __reviseMidiFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3){
	float low=0;
	float high=0;

	int det=0;
	int offset=0;

	if(!isEdge){
		det=2;
		offset=1;
	}

	low=auditory_freToMidi(lowFre)-offset;
	high=low+num-1+det;

	*lowFre3=auditory_midiToFre(low);
	*highFre3=auditory_midiToFre(high);
}

static void __reviseLogFre(int num,float lowFre,float highFre,int binPerOctave,int isEdge,float *lowFre3,float *highFre3){
	float low=0;
	float high=0;

	int det=0;
	int offset=0;

	if(!isEdge){
		det=2;
		offset=1;
	}

	low=auditory_freToLog(lowFre,binPerOctave)-offset;
	high=low+num-1+det;

	*lowFre3=auditory_logToFre(low,binPerOctave);
	*highFre3=auditory_logToFre(high,binPerOctave);
}

static void __reviseLinearFre(int num,float lowFre,float highFre,float detFre,int isEdge,float *lowFre3,float *highFre3){
	float low=0;
	float high=0;

	int det=0;
	int offset=0;

	if(!isEdge){
		det=2;
		offset=1;
	}

	low=roundf(lowFre/detFre)-offset;
	high=low+num-1+det;

	*lowFre3=low*detFre;
	*highFre3=high*detFre;
}

static void __reviseLinspaceFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3){
	float detFre=0;

	if(!isEdge){
		detFre=(highFre-lowFre)/(num-1);

		*lowFre3=lowFre-detFre;
		*highFre3=highFre+detFre;
	}
	else{
		*lowFre3=lowFre;
		*highFre3=highFre;
	}
}

static void __reviseLogspaceFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3){
	float low=0;
	float high=0;

	float det=0;

	if(!isEdge){
		low=auditory_freToLogspace(lowFre);
		high=auditory_freToLogspace(highFre);

		det=(high-low)/(num-1);

		low=low-det;
		high=high+det;

		*lowFre3=auditory_logspaceToFre(low);
		*highFre3=auditory_logspaceToFre(high);
	}
	else{
		*lowFre3=lowFre;
		*highFre3=highFre;
	}
}

// linear scale
float auditory_freToLinear(float fre,float detFre){
	float value=0;

	value=roundf(fre/detFre);
	return value;
}

float auditory_linearToFre(float value,float detFre){
	float fre=0;

	fre=value*detFre;
	return fre;
}

// linspace scale
float auditory_freToLinspace(float fre){

	return fre;
}

float auditory_linspaceToFre(float value){

	return value;
}

// mel scale
// mel=2595*log10(1+fre/700)
float auditory_freToMel(float fre){
	float mel=0;

	mel=2595*log10f(1+fre/700);
	return mel;
}

// fre=700*(10^(mel/2595)-1)
float auditory_melToFre(float mel){
	float fre=0;

	fre=700*(powf(10, mel/2595)-1);
	return fre;
}

// bark scale 巴克
/***
	bark=26.81*fre/(1960+fre)-0.53
	bark=bark+0.15*(2-bark); bark<2
	bark=bark+0.22*(bark-20.1); bark>20.1
****/
float auditory_freToBark(float fre){
	float bark=0;

	bark=26.81*fre/(1960+fre)-0.53;
	if(bark<2){
		bark=bark+0.15*(2-bark);
	}
	else if(bark>20.1){
		bark=bark+0.22*(bark-20.1);
	}

	return bark;
}

/***
	bark=(bark-0.3)/0.85; bark<2
	bark=(bark+4.422)/1.22; bark>20.1
	fre=1960*(bark+0.53)/(26.28-bark);
****/
float auditory_barkToFre(float bark){
	float fre=0;

	if(bark<2){
		bark=(bark-0.3)/0.85;
	}
	else if(bark>20.1){
		bark=(bark+4.422)/1.22;
	}

	fre=1960*(bark+0.53)/(26.28-bark);

	return fre;
}

// erb scale 等效矩形带宽
/***
	erb=A*log10f(1+fre*0.004368)
	A=1000*logf(10)/(24.7*4.37) ~=21.3654
****/
float auditory_freToErb(float fre){
	float erb=0;
	float a=0;

	// a=1000*logf(10)/(24.7*4.37);
	a=21.3654;
	erb=a*log10f(1+fre*0.004368);

	return erb;
}

/***
	fre=(10^(erb/A)-1)/0.004368
	A=1000*logf(10)/(24.7*4.37) ~=21.3654
****/
float auditory_erbToFre(float erb){
	float fre=0;
	float a=0;

	// a=1000*logf(10)/(24.7*4.37);
	a=21.3654;
	fre=(powf(10, erb/a)-1)/0.004368;

	return fre;
}

// midi sacle
/***
	midi=12*log2(fre/440)+69
****/
float auditory_freToMidi(float fre){
	float midi=0;

	midi=roundf(12*log2(fre/440)+69);

	return midi;
}

float auditory_midiToFre(float midi){
	float fre=0;

	fre=powf(2,(midi-69)/12)*440;

	return fre;
}

float auditory_freToLog(float fre,float binPerOctave){
	float value=0;
	float n=0;

	n=binPerOctave/12;
	value=roundf(binPerOctave*log2(fre/440));

	return value;
}

float auditory_logToFre(float value,float binPerOctave){
	float fre=0;
	float n=0;

	n=binPerOctave/12;
	fre=pow(2,value/binPerOctave)*440;

	return fre;
}

// logspace scale
float auditory_freToLogspace(float fre){
	float value=0;

	value=log2(fre/440);
	return value;
}

float auditory_logspaceToFre(float value){
	float fre=0;

	fre=pow(2, value)*440;
	return fre;
}

void auditory_reviseLogFre(int num,float lowFre,float highFre,int binPerOctave,int isEdge,float *lowFre3,float *highFre3){

	__reviseLogFre(num,lowFre,highFre,binPerOctave,isEdge,lowFre3,highFre3);
}

void auditory_reviseLinearFre(int num,float lowFre,float highFre,float detFre,int isEdge,float *lowFre3,float *highFre3){

	__reviseLinearFre(num,lowFre,highFre,detFre,isEdge,lowFre3,highFre3);
}

void auditory_reviseLinspaceFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3){

	__reviseLinspaceFre(num, lowFre, highFre, isEdge, lowFre3, highFre3);
}

void auditory_reviseLogspaceFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3){

	__reviseLogspaceFre(num,lowFre, highFre, isEdge, lowFre3, highFre3);
}






















