// clang 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "chroma_filterBank.h"

// filterBank相关
void chroma_stftFilterBank(int num,int fftLength,int samplate,
						float *octaveCenter,float *octaveWidth,
						float *mFilterBankArr){
	float center=5;
	float width=2;

	float *octArr=NULL;
	float *widthArr=NULL;

	float *mArr1=NULL;
	float *mArr2=NULL;

	float *vArr1=NULL;

	int half=0;
	int n=0;

	float baseFre=440.0; // 523.251131 ???

	if(num%12!=0||num<12){
		printf("num is error\n");
	}

	if(octaveCenter){
		if(*octaveCenter>0){
			center=*octaveCenter;
		}
	}

	if(octaveWidth){
		if(*octaveWidth>0){
			width=*octaveWidth;
		}
	}

	n=num/12;

	octArr=__vnew(fftLength, NULL);
	widthArr=__vnew(fftLength, NULL);

	mArr1=__vnew(num*fftLength, NULL);
	mArr2=__vnew(num*fftLength, NULL);

	vArr1=__vnew(fftLength, NULL);

	// 1. octArr&widthArr
	for(int i=1;i<fftLength;i++){
		float _fre=0;

		_fre=1.0*i/fftLength*samplate;
		octArr[i]=num*logf(_fre/(baseFre/16))/logf(2.0);
	}
	octArr[0]=octArr[1]-1.5*num;

	for(int i=1;i<fftLength;i++){
		float _value=0;

		_value=octArr[i]-octArr[i-1];
		widthArr[i-1]=(_value>1?_value:1);
	}
	widthArr[fftLength-1]=1;

	// 2. chromaMatrix
	half=roundf(num/2.0);
	for(int i=0;i<num;i++){
		for(int j=0;j<fftLength;j++){
			float _value1=0;
			float _value2=0;
			int _k=0;

			_value1=octArr[j]-i;

			_value1=_value1+half+10*num;
			_k=floorf(_value1/num);
			_value2=_value1-_k*num;

			_value2=_value2-half;
			mArr1[i*fftLength+j]=_value2;
		}
	}

	// 3. guass
	for(int i=0;i<num;i++){
		for(int j=0;j<fftLength;j++){
			float _value1=0;
			float _value2=0;

			_value1=mArr1[i*fftLength+j];

			_value1=2*_value1/widthArr[j];
			_value1=-0.5*_value1*_value1;

			_value2=expf(_value1);
			mArr1[i*fftLength+j]=_value2;
			mArr2[i*fftLength+j]=_value2*_value2;
		}
	}

	// 4. norm
	__msum(mArr2, num, fftLength, 0, vArr1);
	for(int i=0;i<num;i++){
		for(int j=0;j<fftLength;j++){
			float _value1=0;
			float _value2=0;

			_value1=mArr1[i*fftLength+j];
			_value2=sqrtf(vArr1[j]);

			_value2=_value1/_value2;
			mArr1[i*fftLength+j]=_value2;
		}
	}

	// 5. scale
	if(width>0){
		for(int i=0;i<num;i++){
			for(int j=0;j<fftLength/2+1;j++){
				float _value1=0;
				float _value2=0;

				_value1=mArr1[i*fftLength+j];

				_value2=(octArr[j]/num-center)/width;
				_value2=-0.5*_value2*_value2;
				_value2=expf(_value2);

				_value2=_value1*_value2;

				mArr1[i*(fftLength/2+1)+j]=_value2;
				// mFilterBankArr[i*(fftLength/2+1)+j]=_value2;
			}
		}
	}

	// 6. offset ???
	for(int i=3*n,k=0;i<num;i++,k++){
		for(int j=0;j<fftLength/2+1;j++){

			mFilterBankArr[k*(fftLength/2+1)+j]=mArr1[i*(fftLength/2+1)+j];
		}
	}

	for(int i=0,k=num-3*n;i<3*n;i++,k++){
		for(int j=0;j<fftLength/2+1;j++){

			mFilterBankArr[k*(fftLength/2+1)+j]=mArr1[i*(fftLength/2+1)+j];
		}
	}

	free(octArr);
	free(widthArr);

	free(mArr1);
	free(mArr2);

	free(vArr1);
}

/***
	num*cqtLength
	binPerOctave map num
	num<=binPerOctave
****/
void chroma_cqtFilterBank(int num,int cqtLength,int binPerOctave,
						float *minFre,
						float *mFilterBankArr){
	int n=0;
	int offset=0;
	int sub=0;

	int start=0;

	float _minFre=32.703196;
	int midiIndex=0;

	float *arr=NULL;

	if(num>binPerOctave||binPerOctave%num!=0){
		printf("num and binPerOctave not map!!!");
		return;
	}

	if(minFre){
		if(*minFre>0){
			_minFre=*minFre;
		}
	}

	n=binPerOctave/num;
	offset=ceilf(n/2.0);
	sub=n-offset;

	midiIndex=roundf(12*log2(_minFre/440)+69);
	midiIndex=midiIndex%12;
	if(midiIndex>6){
		midiIndex=12-midiIndex;
	}

	if(midiIndex){
		arr=__vnew(num*cqtLength, NULL);
	}
	else{
		arr=mFilterBankArr;
	}

	for(int i=0;i<num;i++){
		if(i){
			start=offset+(i-1)*n;
		}
		for(int j=0;j<cqtLength;j++){
			int _mod=0;

			_mod=j%binPerOctave;
			if(i!=0){
				if(_mod>=start&&_mod<start+n){
					arr[i*cqtLength+j]=1;
				}
				
			}
			else{ // first
				if(_mod>=0&&_mod<offset){
					arr[i*cqtLength+j]=1;
				}

				if(sub){
					if(_mod>=binPerOctave-sub&&_mod<binPerOctave){
						arr[i*cqtLength+j]=1;
					}
				}
			}
		}
	}

	if(midiIndex){
		n=num/binPerOctave;
		for(int i=midiIndex*n,k=0;i<num;i++,k++){
			for(int j=0;j<cqtLength;j++){

				mFilterBankArr[k*cqtLength+j]=arr[i*cqtLength+j];
			}
		}

		for(int i=0,k=num-midiIndex*n;i<midiIndex*n;i++,k++){
			for(int j=0;j<cqtLength;j++){

				mFilterBankArr[k*cqtLength+j]=arr[i*cqtLength+j];
			}
		}

		free(arr);
	}
}

void chroma_genericFilterBank(int num,int fftLength,int samplate,
							int freLength,float *freBandArr,
							float *mFilterBankArr){
	
}






















