// 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_complex.h"

#include "flux_spectral.h"

static void __spectral_pd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
						int *indexArr,int indexLength,
						int isWeight,int isNorm,
						float *vArr);

static void __spectral_cd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
						int *indexArr,int indexLength,
						int isRectify,
						float *vArr);

void spectral_flatness(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,
					float *vArr){
	double n1=1;
	float m1=0;

	for(int i=0;i<nLength;i++){
		// n1=1;	
		// for(int j=start;j<=end;j++){
		// 	n1*=mDataArr[i*mLength+j];
		// }

		// convert log
		n1=0;	
		for(int j=0;j<indexLength;j++){
			int _index=0;

			_index=indexArr[j];
			n1+=logf(mDataArr[i*mLength+_index]+2.0e-16);
		}

		// n1=powf(n1, 1.0/(end-start));
		n1=n1/indexLength;
		n1=expf(n1);
		m1=sumArr[i]/indexLength;

		if(m1){
			vArr[i]=n1/m1; // m1 eps ???
		}
		else{
			vArr[i]=0;
		}
	}
}

// step>=1 isExp 0 type 0 sum 1 mean 
void spectral_flux(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					int step,float p,int isPostive,int isExp,int type,
					float *vArr){
	float value=0;

	if(step<1){
		step=1;
	}

	memset(vArr, 0, sizeof(float )*step);
	for(int i=step;i<nLength;i++){
		value=0;
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=mDataArr[i*mLength+_index]-mDataArr[(i-step)*mLength+_index];
			if(isPostive){
				v1=(v1>0?v1:0);
			}
			else{
				v1=fabsf(v1);
			}

			if(p==2.0){
				v1*=v1;
			}
			else{
				v1=powf(v1, p);
			}
			
			value+=v1;
		}

		if(type){ // mean
			value/=indexLength;
		}

		if(isExp){
			vArr[i]=powf(value, 1.0/p); 
		}
		else{
			vArr[i]=value;
		}
	}
}

void spectral_rolloff(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,float threshold,
					float *vArr){
	float n1=0;
	float m1=0;

	int index=0;

	for(int i=0;i<nLength;i++){
		n1=0;	
		m1=sumArr[i]*threshold;
		for(int j=0;j<indexLength;j++){
			float v=0;
			int _index=0;

			_index=indexArr[j];
			v=fabsf(mDataArr[i*mLength+_index]);
			n1+=v;
			if(n1>=m1){
				// float d1=0;
				// float d2=0;

				// d1=n1-m1;
				// d2=m1-(n1-v);
				// if(d1<d2){
				// 	index=j;
				// }
				// else{
				// 	index=(j-1<start?start:j-1);
				// }

				index=_index;
				break;
			}
		}

		vArr[i]=freArr[index]; // m1 eps ???
	}
}

void spectral_centroid(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,
					float *vArr){
	float n1=0;
	float m1=0;

	for(int i=0;i<nLength;i++){
		n1=0;	
		m1=sumArr[i];
		for(int j=0;j<indexLength;j++){
			int _index=0;

			_index=indexArr[j];
			n1+=freArr[_index]*mDataArr[i*mLength+_index];
			// m1+=mDataArr[i*mLength+j];
		}

		if(m1){
			vArr[i]=n1/m1; // m1 eps ???
		}
		else{
			vArr[i]=0;
		}
	}
}

void spectral_spread(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,float *cArr,
					float *vArr){
	float n1=0;
	float m1=0;

	for(int i=0;i<nLength;i++){
		n1=0;	
		m1=sumArr[i];
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=(freArr[_index]-cArr[i]);
			n1+=v1*v1*
				mDataArr[i*mLength+_index];
		}

		if(m1){
			vArr[i]=sqrtf(n1/m1); // m1 eps ???
		}
		else{
			vArr[i]=0;
		}
	}
}

void spectral_skewness(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,float *cArr1,float *cArr2,
					float *vArr){
	float n1=0;
	float m1=0;

	for(int i=0;i<nLength;i++){
		n1=0;	
		m1=cArr2[i]*cArr2[i]*cArr2[i]*sumArr[i];
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=(freArr[_index]-cArr1[i]);
			n1+=v1*v1*v1*
				mDataArr[i*mLength+_index];
		}

		if(m1){
			vArr[i]=n1/m1; // m1 eps ???
		}
		else{
			vArr[i]=0;
		}
	}
}

void spectral_kurtosis(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,float *cArr1,float *cArr2,
					float *vArr){
	float n1=0;
	float m1=0;

	for(int i=0;i<nLength;i++){
		n1=0;	
		m1=cArr2[i]*cArr2[i]*cArr2[i]*cArr2[i]*sumArr[i];
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=(freArr[_index]-cArr1[i]);
			n1+=v1*v1*v1*v1*
				mDataArr[i*mLength+_index];
		}

		if(m1){
			vArr[i]=n1/m1; // m1 eps ???
		}
		else{
			vArr[i]=0;
		}
	}
}

void spectral_entropy(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *sumArr,int isNorm,
					float *vArr){
	float n1=0;
	float m1=0;

	for(int i=0;i<nLength;i++){
		n1=0;	
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=mDataArr[i*mLength+_index]/sumArr[i];
			n1+=v1*log2f(v1+1e-16);
		}

		if(isNorm){ // matlab
			m1=log2f(indexLength);
			if(m1){
				vArr[i]=-n1/m1; // m1 eps ???
			}
			else{
				vArr[i]=0;
			}
		}
		else{ // song
			vArr[i]=-n1;
		}
	}
}

void spectral_crest(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,
					float *vArr){
	float n1=0;
	float m1=0;

	for(int i=0;i<nLength;i++){
		n1=mDataArr[i*mLength+indexArr[0]];	
		m1=sumArr[i]/indexLength;
		for(int j=1;j<indexLength;j++){
			int _index=0;

			_index=indexArr[j];
			if(n1<mDataArr[i*mLength+_index]){
				n1=mDataArr[i*mLength+_index];
			}
		}

		if(m1){
			vArr[i]=n1/m1; // m1 eps ???
		}
		else{
			vArr[i]=0;
		}
	}
}

void spectral_slope(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *meanFreArr,float *meanValueArr,
					float *vArr){
	float n1=0;
	float m1=0;

	for(int i=0;i<nLength;i++){
		n1=0;	
		m1=0;
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=(freArr[_index]-meanFreArr[i]);
			n1+=v1*(mDataArr[i*mLength+_index]-meanValueArr[i]);
			m1+=v1*v1;
		}

		if(m1){
			vArr[i]=n1/m1; // m1 eps ???
		}
		else{
			vArr[i]=0;
		}
	}
}

void spectral_decrease(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *sumArr,
					float *vArr){
	float n1=0;
	float m1=0;

	for(int i=0;i<nLength;i++){
		n1=0;
		m1=sumArr[i]-mDataArr[i*mLength+indexArr[0]];	
		for(int j=1;j<indexLength;j++){
			int _index=0;

			_index=indexArr[j];
			n1+=(mDataArr[i*mLength+_index]-mDataArr[i*mLength+indexArr[0]])/(_index);
		}

		if(m1){
			vArr[i]=n1/m1; // m1 eps ???
		}
		else{
			vArr[i]=0;
		}
	}
}

void spectral_bandWidth(float *mDataArr,int nLength,int mLength,
						int *indexArr,int indexLength,
						float *freArr,float *cArr,float p,
						float *vArr){
	float value1=0;

	for(int i=0;i<nLength;i++){
		value1=0;	
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=(freArr[_index]-cArr[i]);
			if(p==2.0){
				v1*=v1;
			}
			else{
				v1=powf(v1, p);
			}

			value1+=(mDataArr[i*mLength+_index]*v1);
		}

		if(p!=1.0){
			value1=powf(value1, 1.0/p);
		}
		vArr[i]=value1;
	}
}

void spectral_rms(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr){
	float value1=0;

	for(int i=0;i<nLength;i++){
		value1=0;	
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=mDataArr[i*mLength+_index];
			v1*=v1;
			if(_index==0||
				(mLength%2==0&&_index==mLength-1)){
				
				v1*=0.5;
			}

			value1+=v1;
		}

		value1=sqrtf(2*value1/(mLength*mLength));
		vArr[i]=value1;
	}
}

// ['hfc', 'sd', 'sf', 'mkl', 'pd', 'wpd', 'nwpd', 'cd', 'rcd']
// type 0 sum 1 mean
void spectral_hfc(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr){
	float value=0;

	vArr[0]=0;
	for(int i=0;i<nLength;i++){
		value=0;
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=mDataArr[i*mLength+_index]*_index;
			value+=v1;
		}

		vArr[i]=value;
	}
}

// sum(diff); positive=1
void spectral_sd(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				int step,int isPostive,
				float *vArr){
	float value=0;

	if(step<1){
		step=1;
	}
	
	memset(vArr, 0, sizeof(float )*step);
	for(int i=step;i<nLength;i++){
		value=0;
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=mDataArr[i*mLength+_index]-mDataArr[(i-step)*mLength+_index];
			if(isPostive){
				v1=(v1>0?v1:0);
			}
			else{
				v1=fabsf(v1);
			}

			value+=v1;
		}

		vArr[i]=value;
	}
}

// sum(diff^2); positive=1
void spectral_sf(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				int step,int isPostive,
				float *vArr){
	float value=0;

	if(step<1){
		step=1;
	}

	memset(vArr, 0, sizeof(float )*step);
	for(int i=step;i<nLength;i++){
		value=0;
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=mDataArr[i*mLength+_index]-mDataArr[(i-step)*mLength+_index];
			if(isPostive){
				v1=(v1>0?v1:0);
			}
			else{
				v1=fabsf(v1);
			}

			value+=v1*v1;
		}

		vArr[i]=value; 
	}
}

// type 0 sum 1 mean
void spectral_mkl(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				int type,
				float *vArr){
	float value=0;

	vArr[0]=0;
	for(int i=1;i<nLength;i++){
		value=0;
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=mDataArr[i*mLength+_index]/(mDataArr[(i-1)*mLength+_index]+1e-16);
			// log compress
			v1=logf(1+v1);
			value+=v1;
		}

		if(type){
			value/=indexLength;
		}

		vArr[i]=value;
	}
}

static void __spectral_pd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
						int *indexArr,int indexLength,
						int isWeight,int isNorm,
						float *vArr){
	float value=0;

	vArr[0]=0;
	for(int i=2;i<nLength;i++){
		value=0;
		for(int j=0;j<indexLength;j++){
			float v1=0;
			int _index=0;

			_index=indexArr[j];
			v1=mPhaseArr[i*mLength+_index]-2*mPhaseArr[(i-1)*mLength+_index]+mPhaseArr[(i-2)*mLength+_index];
			// 0~2pi --> -pi~pi
			
			v1=fabsf(v1);
			if(isWeight||isNorm){
				v1=v1*mSpecArr[i*mLength+_index];
			}

			value+=v1;
		}

		value/=indexLength;
		if(isNorm){
			float _m1=0;

			for(int j=0;j<indexLength;j++){
				int _index=0;

				_index=indexArr[j];
				_m1+=mSpecArr[i*mLength+_index];
			}
			_m1=_m1/indexLength;

			value=value/(_m1+1e-16);
		}

		vArr[i]=value;
	}
}

void spectral_pd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr){
	
	__spectral_pd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				0,0,
				vArr);
}

void spectral_wpd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr){

	__spectral_pd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				1,0,
				vArr);
}

void spectral_nwpd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr){

	__spectral_pd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				0,1,
				vArr);
}

static void __spectral_cd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
						int *indexArr,int indexLength,
						int isRectify,
						float *vArr){
	float value=0;

	vArr[0]=0;
	for(int i=1;i<nLength;i++){
		value=0;
		for(int j=0;j<indexLength;j++){
			float v1=0;
			float v2=0;

			float real1=0;
			float image1=0;

			float real2=0;
			float image2=0;

			int _index=0;

			_index=indexArr[j];

			if(isRectify){
				if(mSpecArr[i*mLength+_index]<=mSpecArr[(i-1)*mLength+_index]){
					continue;
				}
			}

			v1=mPhaseArr[i*mLength+_index];
			real1=mSpecArr[i*mLength+_index]*cosf(v1);
			image1=mSpecArr[i*mLength+_index]*sin(v1);

			if(i>1){
				v2=2*mPhaseArr[(i-1)*mLength+_index]-mPhaseArr[(i-2)*mLength+_index];
				
				real2=mSpecArr[(i-1)*mLength+_index]*cosf(v2);
				image2=mSpecArr[(i-1)*mLength+_index]*sin(v2);

				real1-=real2;
				image1-=image2;
			}
			
			v1=sqrtf(real1*real1+image1*image1);
			value+=v1;
		}

		vArr[i]=value;
	}
}

void spectral_cd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr){

	__spectral_cd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				0,
				vArr);
}

void spectral_rcd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr){

	__spectral_cd(mSpecArr,mPhaseArr,nLength,mLength,
				indexArr,indexLength,
				1,
				vArr);
}

// threshold 0/3/6...
void spectral_broadband(float *mDataArr,int nLength,int mLength,
						int *indexArr,int indexLength,
						float threshold,
						float *vArr){
	vArr[0]=0;
	for(int i=1;i<nLength;i++){
		for(int j=0;j<indexLength;j++){
			float diff=0;
			int _index=0;

			_index=indexArr[j];
			diff=10.0*log10f(mDataArr[i*mLength+_index]/(mDataArr[(i-1)*mLength+_index]));
			if(diff>threshold){
				vArr[i]++;
			}
		}
	}
}

/***
	step >=1
	threshold 0
	methodType 'sub'
	dataType 'value'
****/
void spectral_novelty(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					int step,float threshold,
					SpectralNoveltyMethodType *methodType,SpectralNoveltyDataType *dataType,
					float *vArr){
	float value=0;

	SpectralNoveltyMethodType mType=SpectralNoveltyMethod_Sub;
	SpectralNoveltyDataType dType=SpectralNoveltyData_Value;

	if(methodType){
		mType=*methodType;
	}

	if(dataType){
		dType=*dataType;
	}

	if(step<1){
		step=1;
	}

	memset(vArr, 0, sizeof(float )*step);
	for(int i=step;i<nLength;i++){
		value=0;
		for(int j=0;j<indexLength;j++){
			float pre=0;
			float cur=0;

			float v1=0;
			int _index=0;

			_index=indexArr[j];

			cur=mDataArr[i*mLength+_index];
			pre=mDataArr[(i-step)*mLength+_index];

			if(mType==SpectralNoveltyMethod_Sub){
				v1=cur-pre;
			}
			else if(mType==SpectralNoveltyMethod_Entroy){
				v1=logf(cur/(pre+1e-16));
			}
			else if(mType==SpectralNoveltyMethod_KL){
				v1=cur*logf(cur/(pre+1e-16));
			}
			else{ // IS
				v1=cur/(pre+1e-16)-logf(cur/(pre+1e-16))-1;
			}

			if(dType==SpectralNoveltyData_Value){
				if(v1>threshold){
					value+=v1;
				}
			}
			else{ // Number
				if(v1>threshold){
					value++;
				}
			}
		}

		vArr[i]=value;
	}
}

void spectral_energy(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					int isPower,int isLog,float gamma,
					float *vArr){
	for(int i=0;i<nLength;i++){
		vArr[i]=0;
		for(int j=0;j<indexLength;j++){
			float _value=0;
			int _index=0;

			_index=indexArr[j];
			_value=mDataArr[i*mLength+_index];
			if(!isPower){
				_value*=_value;
			}

			if(isLog){
				if(gamma<=0){
					gamma=10;
				}

				_value=logf(1+gamma*_value);
			}

			vArr[i]+=_value;
		}

		vArr[i]/=indexLength; // fre domain ???
	}	
}










