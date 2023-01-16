// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "../dsp/flux_window.h"
#include "../dsp/fft_algorithm.h"
#include "../dsp/dft_algorithm.h"

#include "auditory_filterBank.h"
#include "cwt_filterBank.h"

static void __cwt_morseFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,float cf,
								float *mFilterBankArr);

static void __cwt_morletFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr);

static void __cwt_bumpFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr);

static void __cwt_paulFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr);

static void __cwt_dogFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr);

static void __cwt_mexicanFilterBank(float *wArr,int wLength,float *sArr,int sLength,
									float gamma,float beta,
									float *mFilterBankArr);

static void __cwt_poissonFilterBank(float *wArr,int wLength,float *sArr,int sLength,
									float gamma,float beta,
									float *mFilterBankArr);

static void __cwt_hermitFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr);

static void __cwt_rickerFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr);

static void __cwt_calBandEdge(int num,int dataLength,int samplate,
							float lowFre,float highFre,int isEdge,
							void *callback1,void *callback2,float ref,
							float **freBandArr,int **binBandArr);

static float __calMorseCenterFre(float gamma,float beta);
static float __calMorseSD(float gamma,float beta);

static long double __calPaulFactor(float m);
static long double __calDogFactor(float m);
static long double __calPoissonFactor(float m);
static long double __calHermitFactor(float gamma);

/***
	cwt-filterBank相关 

	log+linear+mel/bark/erb 都不包含边界
		log 针对cwt-spectrogram形式上包含
		log is use default
		
	Morse gamma=3 beta=20; 
	Morlet gamma=6 beta=2; 
	Bump gamma=5 beta=0.6;
	
	Paul gamma=4
	DOG gamma=2 beta=2;
	Mexih beta=2

	type semantic is wavelet type and window type(styleType)
****/
void cwt_filterBank(int num,int dataLength,int samplate,int padLength,
					WaveletContinueType type,float gamma,float beta,
					SpectralFilterBankScaleType scaleType,
					float lowFre,float highFre,int binPerOctave,
					float *mFilterBankArr,
					float *freBandArr,
					int *binBandArr){
	float *wArr=NULL;
	float *sArr=NULL;

	int wLength=0;
	int sLength=0;

	float cf=0;

	int isEdge=0;
	int offset=0;

	float *fArr=NULL;
	int *bArr=NULL;

	// UniFunc||UniFunc1
	void *func1=NULL;
	void *func2=NULL;

	float ref=0;

	isEdge=0;
	if(!isEdge){
		offset=1;
	}

	// isPad in cwt
	wLength=dataLength+2*padLength;
	sLength=num;

	// 0. revise log/linear!!!
	if(scaleType==SpectralFilterBankScale_Octave){ 
		if(binPerOctave&&binPerOctave>=4&&binPerOctave<=48){
			ref=binPerOctave;
		}
		else{
			ref=12;
		}

		auditory_reviseLogFre(num,lowFre,highFre,ref,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Linear){
		ref=samplate*1.0/dataLength;

		auditory_reviseLinearFre(num,lowFre,highFre,ref,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Linspace){
		auditory_reviseLinspaceFre(num,lowFre,highFre,isEdge,&lowFre,&highFre);
	}
	else if(scaleType==SpectralFilterBankScale_Log){
		auditory_reviseLogspaceFre(num,lowFre,highFre,isEdge,&lowFre,&highFre);
	}

	if(type==WaveletContinue_Morse){
		cf=__calMorseCenterFre(gamma,beta);
	}
	else if(type==WaveletContinue_Morlet||
			type==WaveletContinue_Bump){
		cf=gamma;
	}
	else if(type==WaveletContinue_Paul){
		cf=gamma+0.5;
	}
	else if(type==WaveletContinue_DOG){
		cf=sqrtf(gamma+0.5);
	}
	else if(type==WaveletContinue_Mexican){
		cf=sqrtf(2+0.5);
	}
	else if(type==WaveletContinue_Hermit){
		cf=gamma+1;
	}
	else if(type==WaveletContinue_Ricker){
		cf=gamma;
	}
	
	// else if(type==WaveletContinue_Morlet){
	// 	cf=6.0;
	// }
	// else if(type==WaveletContinue_Bump){
	// 	cf=5.0;
	// }

	// 1. cal freBand/binBand
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
	
	__cwt_calBandEdge(num,dataLength,samplate,
					lowFre,highFre,isEdge,
					func1,func2,ref,
					&fArr,&bArr);

	// 2. cal wArr,sArr
	wArr=__vnew(wLength, NULL);
	for(int i=0;i<=wLength/2;i++){
		wArr[i]=i*2*M_PI/wLength;
	}
	for(int i=wLength/2+1,j=wLength/2-1;i<wLength&&j>=0;i++,j--){
		wArr[i]=-wArr[j];
	}

	sArr=__vnew(num, NULL);
	for(int i=sLength-1+offset,j=0;i>=offset;i--,j++){
		float _value=0;

		_value=fArr[i];
		if(_value<1e-6){
			_value=1e-6;
		}
		
		sArr[j]=cf/(_value/samplate*2*M_PI);
	}

	// 3. filterBankArr
	if(type==WaveletContinue_Morse){
		__cwt_morseFilterBank(wArr,wLength,sArr,sLength,
							gamma,beta,cf,
							mFilterBankArr);
	}
	else if(type==WaveletContinue_Morlet){
		__cwt_morletFilterBank(wArr,wLength,sArr,sLength,
							gamma,beta,
							mFilterBankArr);
	}
	else if(type==WaveletContinue_Bump){
		__cwt_bumpFilterBank(wArr,wLength,sArr,sLength,
							gamma,beta,
							mFilterBankArr);
	}
	else if(type==WaveletContinue_Paul){
		__cwt_paulFilterBank(wArr,wLength,sArr,sLength,
							gamma,beta,
							mFilterBankArr);
	}
	else if(type==WaveletContinue_DOG){
		__cwt_dogFilterBank(wArr,wLength,sArr,sLength,
							gamma,beta,
							mFilterBankArr);
	}
	else if(type==WaveletContinue_Mexican){
		__cwt_mexicanFilterBank(wArr,wLength,sArr,sLength,
							gamma,beta,
							mFilterBankArr);
	}
	else if(type==WaveletContinue_Hermit){
		__cwt_hermitFilterBank(wArr,wLength,sArr,sLength,
							gamma,beta,
							mFilterBankArr);
	}
	else if(type==WaveletContinue_Ricker){
		__cwt_rickerFilterBank(wArr,wLength,sArr,sLength,
							gamma,beta,
							mFilterBankArr);
	}

	if(freBandArr){
		memcpy(freBandArr, fArr+offset, sizeof(float )*num);
	}
	
	if(binBandArr){
		memcpy(binBandArr, bArr+offset, sizeof(int )*num);
	}

	// printf("mFilterBankArr is \n");
	// __mdebug(mFilterBankArr, sLength, wLength, 1);
	// printf("\n");

	// printf("freArr is \n");
	// __vdebug(fArr+1, num, 1);
	// printf("\n");

	free(fArr);
	free(bArr);
}

static void __cwt_calBandEdge(int num,int dataLength,int samplate,
							float lowFre,float highFre,int isEdge,
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

	if(!isEdge){ // ! linear
		det=2;
	}

	if(callback1&&callback2){
		if(!ref){
			func1=(UniFunc )callback1;
			func2=(UniFunc )callback2;

			low=func1(lowFre);
			high=func1(highFre);
		}
		else{
			func3=(UniFunc1 )callback1;
			func4=(UniFunc1 )callback2;

			low=func3(lowFre,ref);
			high=func3(highFre,ref);
		}
	}
	else{
		low=lowFre;
		high=highFre;
	}

	// 1. mel/barkArr
	fArr=__vlinspace(low, high, num+det, 0);

	// 2. mel/barkArr=>freArr
	if(callback1&&callback2){
		if(!ref){
			__vmap(fArr, num+det, func2, NULL);
		}
		else{
			__vmap1(fArr, num+det, func4,ref, NULL);
		}
	}
	
	if(freBandArr){
		*freBandArr=fArr;
	}

	// 3. freArr=>binArr vector???
	if(binBandArr){
		bArr=__vnewi(num+det, NULL);
		for(int i=0;i<num+det;i++){
			bArr[i]=roundf(dataLength*fArr[i]/samplate);
		}

		*binBandArr=bArr;
	}

	if(!freBandArr){
		free(fArr);
	}

}

static void __cwt_morseFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,float cf,
								float *mFilterBankArr){
	float factor=0;
	int totaLength=0;

	totaLength=sLength*wLength;
	factor=expf(-beta*logf(cf)+powf(cf, gamma));

	__mdot(sArr, wArr, sLength, 1, 1, wLength, mFilterBankArr);
	for(int i=0;i<totaLength;i++){
		float value1=0;
		float value2=0;

		if(mFilterBankArr[i]>0){
			value1=fabsf(mFilterBankArr[i]);
			if(gamma==3){
				value2=value1*value1*value1;
			}
			else{
				value2=powf(value1, gamma);
			}

			mFilterBankArr[i]=2*factor*expf(beta*logf(value1)-value2);
		}
		else{
			mFilterBankArr[i]=0;
		}
	}
}

static void __cwt_morletFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr){
	int totaLength=0;
	float cf=6;

	cf=gamma;
	totaLength=sLength*wLength;

	__mdot(sArr, wArr, sLength, 1, 1, wLength, mFilterBankArr);
	for(int i=0;i<totaLength;i++){
		float value1=0;

		if(mFilterBankArr[i]>0){
			value1=mFilterBankArr[i];

			value1=-(value1-cf)*(value1-cf)/beta;
			if(mFilterBankArr[i]<=0){
				value1=0;
			}

			mFilterBankArr[i]=2*expf(value1);
		}
		else{
			mFilterBankArr[i]=0;
		}
	}
}

static void __cwt_bumpFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr){
	int totaLength=0;
	
	float cf=5;
	float sigma=0.6;
	float eps=1e-6;

	cf=gamma;
	sigma=beta;
	totaLength=sLength*wLength;

	__mdot(sArr, wArr, sLength, 1, 1, wLength, mFilterBankArr);
	for(int i=0;i<totaLength;i++){
		float value1=0;
		float value2=0;
		float value3=0;

		value1=(mFilterBankArr[i]-cf)/sigma;

		value2=value1*value1;
		value2=-1/(1-value2);

		if(fabsf(value1)<1-eps){
			value3=2*M_E*expf(value2);
			if(isnan(value3)){
				mFilterBankArr[i]=0;
			}
			else{
				mFilterBankArr[i]=value3;
			}
		}
		else{
			mFilterBankArr[i]=0;
		}
	}
}

static void __cwt_paulFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr){
	long double factor=0;
	int totaLength=0;

	totaLength=sLength*wLength;
	factor=__calPaulFactor(gamma);

	__mdot(sArr, wArr, sLength, 1, 1, wLength, mFilterBankArr);
	for(int i=0;i<totaLength;i++){
		long double value1=0;

		if(mFilterBankArr[i]>0){
			value1=mFilterBankArr[i];

			value1=powl(value1, gamma)*expl(-value1);
			mFilterBankArr[i]=factor*value1;
		}
		else{
			mFilterBankArr[i]=0;
		}
	}
}

static void __cwt_dogFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr){
	long double factor=0;
	int totaLength=0;

	totaLength=sLength*wLength;
	factor=__calDogFactor(gamma);

	__mdot(sArr, wArr, sLength, 1, 1, wLength, mFilterBankArr);
	for(int i=0;i<totaLength;i++){
		long double value1=0;

		if(mFilterBankArr[i]>0){
			value1=mFilterBankArr[i];
			value1=powl(value1, gamma)*expl(-value1*value1/beta);

			mFilterBankArr[i]=factor*value1;
		}
		else{
			mFilterBankArr[i]=0;
		}
	}
}

static void __cwt_mexicanFilterBank(float *wArr,int wLength,float *sArr,int sLength,
									float gamma,float beta,
									float *mFilterBankArr){
	
	__cwt_dogFilterBank(wArr,wLength,sArr,sLength,
						2,beta,
						mFilterBankArr);
}

static void __cwt_poissonFilterBank(float *wArr,int wLength,float *sArr,int sLength,
									float gamma,float beta,
									float *mFilterBankArr){
	long double factor=0;
	int totaLength=0;

	totaLength=sLength*wLength;
	factor=__calPoissonFactor(gamma);

	__mdot(sArr, wArr, sLength, 1, 1, wLength, mFilterBankArr);
	for(int i=0;i<totaLength;i++){
		long double value1=0;

		if(mFilterBankArr[i]>0){
			value1=mFilterBankArr[i];

			value1=powl(value1, gamma)*expl(-value1);
			mFilterBankArr[i]=factor*value1;
		}
		else{
			mFilterBankArr[i]=0;
		}
	}
}

static void __cwt_hermitFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr){
	long double factor=0;
	int totaLength=0;

	totaLength=sLength*wLength;
	factor=__calHermitFactor(gamma);

	__mdot(sArr, wArr, sLength, 1, 1, wLength, mFilterBankArr);
	for(int i=0;i<totaLength;i++){
		long double value1=0;

		if(mFilterBankArr[i]>0){
			value1=mFilterBankArr[i];

			value1=(value1-gamma)*(1+value1-gamma)*expl(-(value1-gamma)*(value1-gamma)/beta);
			mFilterBankArr[i]=factor*value1;
		}
		else{
			mFilterBankArr[i]=0;
		}
	}
}

static void __cwt_rickerFilterBank(float *wArr,int wLength,float *sArr,int sLength,
								float gamma,float beta,
								float *mFilterBankArr){

	long double factor=0;
	int totaLength=0;

	totaLength=sLength*wLength;
	factor=2.0/sqrtf(M_PI);

	__mdot(sArr, wArr, sLength, 1, 1, wLength, mFilterBankArr);
	for(int i=0;i<totaLength;i++){
		long double value1=0;

		if(mFilterBankArr[i]>0){
			value1=mFilterBankArr[i];

			value1=value1*value1/(gamma*gamma*gamma)*expl(-value1*value1/(gamma*gamma));
			mFilterBankArr[i]=factor*value1;
		}
		else{
			mFilterBankArr[i]=0;
		}
	}
}

// 2/sqrt(gamma)pi^(-0.25)
static long double __calHermitFactor(float gamma){
	long double factor=0;

	
	factor=2.0/sqrtf(gamma)*powl(M_PI, -0.25);

	return factor;
}

// 1/gamma(m+1)
static long double __calPoissonFactor(float m){
	long double factor=0;

	
	factor=1.0/util_gamma(m+1, NULL);

	return factor;
}

// 2^m/sqrt(m*prod(2:(2*m-1))))
static long double __calPaulFactor(float m){
	long double factor=0;

	long double value1=0;
	int p=0;

	value1=1;
	p=roundf(m); 
	for(int i=2*p-1;i>=2;i--){
		value1*=i;
	}

	factor=powl(2, p)/sqrtl(p*value1);

	return factor;
}

// -((1i^m)/sqrt(gamma(m+0.5)))
static long double __calDogFactor(float m){
	long double factor=0;

	long double value1=0;
	int p=0;

	p=roundf(m); 
	value1=-1.0/sqrtl(util_gammal(p+0.5, NULL));

	factor=value1;
	if((p/2)%2==1){
		factor=-factor;
	}

	return factor;
}

// beta/(gamma)^(1/gamma)
static float __calMorseCenterFre(float gamma,float beta){
	float cf=0;

	cf=expf(1.0/gamma*(logf(beta)-logf(gamma)));

	return cf;
}

static float __logmorse(float gamma,float beta){

	return beta/gamma*(1+logf(gamma)-logf(beta));
}

static float __calMorseSD(float gamma,float beta){
	float sd=0;

	float ra=0;
	float rb=0;
	float rc=0;

	float loga=0;
	float logb=0;
	float logc=0;

	float v1=0;
	float v2=0;

	v1=__logmorse(gamma, beta);
	v2=__logmorse(gamma, 2*beta);

	ra=2*v1-
		2*__logmorse(gamma, beta-1)+
		__logmorse(gamma, 2*(beta-1))-
		v2;

	rb=2*v1-
		2*__logmorse(gamma, beta-1+gamma)+
		__logmorse(gamma, 2*(beta-1+gamma))-
		v2;

	rc=2*v1-
		2*__logmorse(gamma, beta-1+gamma/2)+
		__logmorse(gamma, 2*(beta-1+gamma/2))-
		v2;

	loga=ra+
			2/gamma*logf(beta/gamma)+
			2*logf(beta)+
			logf(util_gamma((2*(beta-1)+1)/gamma, NULL))-
			logf(util_gamma((2*beta+1)/gamma,NULL));

	logb=rb+
			2/gamma*logf(beta/gamma)+
			2*logf(gamma)+
			logf(util_gamma((2*(beta-1+gamma)+1)/gamma, NULL))-
			logf(util_gamma((2*beta+1)/gamma,NULL));

	logc=rc+
			2/gamma*logf(beta/gamma)+
			logf(2)+logf(beta)+logf(gamma)+
			logf(util_gamma((2*(beta-1+gamma/2)+1)/gamma, NULL))-
			logf(util_gamma((2*beta+1)/gamma,NULL));	

	sd=sqrtf(expf(loga)+expf(logb)-expf(logc));	

	return sd;
}







