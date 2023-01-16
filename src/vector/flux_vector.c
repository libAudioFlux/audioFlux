// 

#include <string.h>
#include <math.h>

#include "flux_vector.h"

static float __arr_max(float *vArr,int length);
static void __mREZ(float *mArr1,int nLength,int mLength,int axis,float *vArr3,float (*func)(float *,int ));
static void __vunwrap1(float **vArr1,int length);

// 针对 type 0即符合乘法
void __mdot(float *mArr1,float *mArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			float *mArr3){
	if(mLength1!=nLength2){
		return;
	}

	for(int i=0;i<nLength1;i++){
		for(int j=0;j<mLength2;j++){
			double _value=0;

			for(int k=0;k<mLength1;k++){ // arr1.row*arr2.col
				_value+=mArr1[i*mLength1+k]*mArr2[j+k*mLength2];
			}

			mArr3[i*mLength2+j]=_value;
		}
	}
}

// 针对 type 1即m2转置符合乘法 
void __mdot1(float *mArr1,float *mArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			float *mArr3){
	if(mLength1!=mLength2){
		return;
	}

	for(int i=0;i<nLength1;i++){
		for(int j=0;j<nLength2;j++){
			double _value=0;

			for(int k=0;k<mLength1;k++){ // arr1.row*arr2.row
				_value+=mArr1[i*mLength1+k]*mArr2[j*mLength2+k];
			}

			mArr3[i*nLength2+j]=_value;
		}
	}
}

/***
	n1*m1 n2*m2 => n3*m3
	type
		0 n1*m2 m1==n2
		1 n1*n2 m1==m2
		2 m1*m2 n1==n2
		3 m1*n2 n1==m2
	不能in place操作
	return 0 sucess 1 error
****/
int __mdot2(float *mArr1,float *mArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			int *type,
			float *mArr3){
	int status=0;

	int nLen3=0;
	int mLen3=0;

	int _type=0;

	if(type){
		_type=*type;
	}

	if(_type==0){ // n1*m2
		nLen3=nLength1;
		mLen3=mLength2;
		if(mLength1!=nLength2){
			return -1;
		}
	}
	else if(_type==1){ // n1*n2
		nLen3=nLength1;
		mLen3=nLength2;
		if(mLength1!=mLength2){
			return -1;
		}
	}
	else if(_type==2){ // m1*m2
		nLen3=mLength1;
		mLen3=mLength2;
		if(nLength1!=nLength2){
			return -1;
		}
	}
	else{ // m1*n2
		nLen3=mLength1;
		mLen3=nLength2;
		if(nLength1!=mLength2){
			return -1;
		}
	}

	for(int i=0;i<nLen3;i++){
		for(int j=0;j<mLen3;j++){
			double _value=0;

			if(_type==0){ // n1*m2
				for(int k=0;k<mLength1;k++){ // arr1.row*arr2.col
					_value+=mArr1[i*mLength1+k]*mArr2[j+k*mLength2];
				}
			}
			else if(_type==1){ // n1*n2 
				for(int k=0;k<mLength1;k++){ // arr1.row*arr2.T.col(arr2.row)
					_value+=mArr1[i*mLength1+k]*mArr2[j*mLength2+k];
				}
			}
			else if(_type==2){ // m1*m2 
				for(int k=0;k<nLength1;k++){ // arr1.T.row(arr1.col)*arr2.col
					_value+=mArr1[i+k*mLength1]*mArr2[j+k*mLength2];
				}
			}
			else{ // m1*n2
				for(int k=0;k<nLength1;k++){ // arr1.T.row(arr1.col)*arr2.T.col(arr2.row)
					_value+=mArr1[i+k*mLength1]*mArr2[j*mLength2+k];
				}
			}

			mArr3[i*mLen3+j]=_value;
		}
	}

	return status;
}

void __msub(float *mArr1,float *mArr2,int nLength,int mLength,float *mArr3){

	for(int i=0;i<nLength*mLength;i++){
		mArr3[i]=mArr1[i]-mArr2[i];
	}
}

/***
	A n*m
	1*m/n*1=>0/1/ => axis 0/1
****/
void __mmul_vector(float *mArr1,float *vArr,int type,int nLength,int mLength,int axis,float *mArr3){
	float *arr=NULL;

	int nLen=0;
	int mLen=0;

	if(mArr3){
		arr=mArr3;
	}
	else{
		arr=mArr1;
	}

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	for(int i=0;i<nLen;i++){
		for(int j=0;j<mLen;j++){
			if(axis==0){ // 1*m
				arr[i+j*mLength]=mArr1[i+j*mLength]*(type==0?vArr[j]:vArr[i]);
			}
			else{ // n*1
				arr[i*mLength+j]=mArr1[i*mLength+j]*(type?vArr[j]:vArr[i]);
			}
		}
	}
}

void __mdiv_vector(float *mArr1,float *vArr,int type,int nLength,int mLength,int axis,float *mArr3){
	float *arr=NULL;

	int nLen=0;
	int mLen=0;

	if(mArr3){
		arr=mArr3;
	}
	else{
		arr=mArr1;
	}

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	for(int i=0;i<nLen;i++){
		for(int j=0;j<mLen;j++){
			if(axis==0){ // 1*m
				if(mArr1[i+j*mLength]){
					arr[i+j*mLength]=mArr1[i+j*mLength]/(type==0?vArr[j]:vArr[i]);
				}
				else{
					arr[i+j*mLength]=0;
				}
			}
			else{ // n*1
				if(mArr1[i*mLength+j]){
					arr[i*mLength+j]=mArr1[i*mLength+j]/(type?vArr[j]:vArr[i]);
				}
				else{
					arr[i*mLength+j]=0;
				}
			}
		}
	}
}

void __madd_vector(float *mArr1,float *vArr,int type,int nLength,int mLength,int axis,float *mArr3){
	float *arr=NULL;

	int nLen=0;
	int mLen=0;

	if(mArr3){
		arr=mArr3;
	}
	else{
		arr=mArr1;
	}

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	for(int i=0;i<nLen;i++){
		for(int j=0;j<mLen;j++){
			if(axis==0){ // 1*m
				arr[i+j*mLength]=mArr1[i+j*mLength]+(type==0?vArr[j]:vArr[i]);
			}
			else{ // n*1
				arr[i*mLength+j]=mArr1[i*mLength+j]+(type?vArr[j]:vArr[i]);
			}
		}
	}
}

void __msub_vector(float *mArr1,float *vArr,int type,int nLength,int mLength,int axis,float *mArr3){
	float *arr=NULL;

	int nLen=0;
	int mLen=0;

	if(mArr3){
		arr=mArr3;
	}
	else{
		arr=mArr1;
	}

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	for(int i=0;i<nLen;i++){
		for(int j=0;j<mLen;j++){
			if(axis==0){ // 1*m
				arr[i+j*mLength]=mArr1[i+j*mLength]-(type==0?vArr[j]:vArr[i]);
			}
			else{ // n*1
				arr[i*mLength+j]=mArr1[i*mLength+j]-(type?vArr[j]:vArr[i]);
			}
		}
	}
}

void __mmul_value(float *mArr1,float value,int nLength,int mLength,float *mArr3){
	float *mArr=NULL;

	if(mArr3){
		mArr=mArr3;
	}
	else{
		mArr=mArr1;
	}

	for(int i=0;i<nLength;i++){
		for(int j=0;j<mLength;j++){
			mArr[i*mLength+j]=mArr1[i*mLength+j]*value;
		}
	}
}

/***
	axis 0 row 1 col -1/other all
****/
void __msum(float *mArr1,int nLength,int mLength,int axis,float *vArr3){
	int nLen=0;
	int mLen=0;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			double _value=0;

			for(int j=0;j<mLen;j++){
				if(axis==0){ // col遍历
					_value+=mArr1[i+j*mLength];
				}
				else{ // row 遍历
					_value+=mArr1[i*mLength+j];
				}
			}

			vArr3[i]=_value;
		}
	}
	else{ 
		vArr3[0]=__vsum(mArr1, nLength*mLength);
	}
}

void __mmin(float *mArr1,int nLength,int mLength,int axis,float *vArr3,int *indexArr3){
	float min=0;
	int index=0;

	int nLen=0;
	int mLen=0;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			for(int j=0;j<mLen;j++){
				if(axis==0){ // col遍历
					if(j==0){
						vArr3[i]=mArr1[i+j*mLength];
						if(indexArr3){
							indexArr3[i]=j;	
						}
					}
					else{
						if(vArr3[i]>mArr1[i+j*mLength]){
							vArr3[i]=mArr1[i+j*mLength];
							if(indexArr3){
								indexArr3[i]=j;
							}
						}
					}
					
				}
				else{ // row遍历
					if(j==0){
						vArr3[i]=mArr1[i*mLength+j];
						if(indexArr3){
							indexArr3[i]=j;
						}
					}
					else{
						if(vArr3[i]>mArr1[i*mLength+j]){
							vArr3[i]=mArr1[i*mLength+j];
							if(indexArr3){
								indexArr3[i]=j;
							}
						}
					}
				}
			}
		}
	}
	else{ 
		index=__vmin(mArr1, nLength*mLength,&min);
		vArr3[0]=min;
		if(indexArr3){
			indexArr3[0]=index;
		}
	}
}

void __mmax(float *mArr1,int nLength,int mLength,int axis,float *vArr3,int *indexArr3){
	float max=0;
	int index=0;

	int nLen=0;
	int mLen=0;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			for(int j=0;j<mLen;j++){
				if(axis==0){ // col遍历
					if(j==0){
						vArr3[i]=mArr1[i+j*mLength];
						if(indexArr3){
							indexArr3[i]=j;	
						}
					}
					else{
						if(vArr3[i]<mArr1[i+j*mLength]){
							vArr3[i]=mArr1[i+j*mLength];
							if(indexArr3){
								indexArr3[i]=j;
							}
						}
					}
					
				}
				else{ // row遍历
					if(j==0){
						vArr3[i]=mArr1[i*mLength+j];
						if(indexArr3){
							indexArr3[i]=j;
						}
					}
					else{
						if(vArr3[i]<mArr1[i*mLength+j]){
							vArr3[i]=mArr1[i*mLength+j];
							if(indexArr3){
								indexArr3[i]=j;
							}
						}
					}
				}
			}
		}
	}
	else{ 
		index=__vmax(mArr1, nLength*mLength,&max);
		vArr3[0]=max;
		if(indexArr3){
			indexArr3[0]=index;
		}
	}
}

void __mmean(float *mArr1,int nLength,int mLength,int axis,float *vArr3){

	if(!vArr3){
		return;
	}

	__msum(mArr1,nLength,mLength,axis,vArr3);
	if(axis==0||axis==1){
		if(axis==0){
			__vdiv_value(vArr3, nLength, mLength, vArr3);
		}
		else{ // axis==1
			__vdiv_value(vArr3, mLength, nLength, vArr3);
		}
	}
	else{
		vArr3[0]/=(nLength*mLength);
	}
}

// 转__v系列
void __mmedian(float *mArr1,int nLength,int mLength,int axis,float *vArr3){
	int nLen=0;
	int mLen=0;

	float *vArr=NULL;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;

			vArr=__vnew(mLen, NULL);
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			if(axis==0){
				for(int j=0;j<mLen;j++){
					if(axis==0){ // col遍历
						vArr[j]=mArr1[i+j*mLength];
					}
				}
			}
			else{
				vArr=mArr1+i*mLength;
			}

			vArr3[i]=__vmedian(vArr, mLen);
		}

		if(axis==0){
			free(vArr);
		}
	}
}

void __mvar(float *mArr1,int nLength,int mLength,int axis,int type,float *vArr3){
	int nLen=0;
	int mLen=0;

	float *vArr=NULL;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;

			vArr=__vnew(mLen, NULL);
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			if(axis==0){
				for(int j=0;j<mLen;j++){
					if(axis==0){ // col遍历
						vArr[j]=mArr1[i+j*mLength];
					}
				}
			}
			else{
				vArr=mArr1+i*mLength;
			}

			vArr3[i]=__vvar(vArr, mLen, type);
		}

		if(axis==0){
			free(vArr);
		}
	}
}

void __mstd(float *mArr1,int nLength,int mLength,int axis,int type,float *vArr3){
	int nLen=0;
	int mLen=0;

	float *vArr=NULL;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;

			vArr=__vnew(mLen, NULL);
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			if(axis==0){
				for(int j=0;j<mLen;j++){
					if(axis==0){ // col遍历
						vArr[j]=mArr1[i+j*mLength];
					}
				}
			}
			else{
				vArr=mArr1+i*mLength;
			}

			vArr3[i]=__vstd(vArr, mLen, type);
		}

		if(axis==0){
			free(vArr);
		}
	}
}

void __mcov(float *mArr1,float *mArr2,int nLength,int mLength,int axis,int type,float *vArr3){
	int nLen=0;
	int mLen=0;

	float *vArr1=NULL;
	float *vArr2=NULL;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;

			vArr1=__vnew(mLen, NULL);
			vArr2=__vnew(mLen, NULL);
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			if(axis==0){
				for(int j=0;j<mLen;j++){
					if(axis==0){ // col遍历
						vArr1[j]=mArr1[i+j*mLength];
						vArr2[j]=mArr2[i+j*mLength];
					}
				}
			}
			else{
				vArr1=mArr1+i*mLength;
				vArr2=mArr2+i*mLength;
			}

			vArr3[i]=__vcov(vArr1,vArr2, mLen, type);
		}

		if(axis==0){
			free(vArr1);
			free(vArr2);
		}
	}
}

void __mcorrcoef(float *mArr1,float *mArr2,int nLength,int mLength,int axis,float *vArr3){
	int nLen=0;
	int mLen=0;

	float *vArr1=NULL;
	float *vArr2=NULL;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;

			vArr1=__vnew(mLen, NULL);
			vArr2=__vnew(mLen, NULL);
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			if(axis==0){
				for(int j=0;j<mLen;j++){
					if(axis==0){ // col遍历
						vArr1[j]=mArr1[i+j*mLength];
						vArr2[j]=mArr2[i+j*mLength];
					}
				}
			}
			else{
				vArr1=mArr1+i*mLength;
				vArr2=mArr2+i*mLength;
			}

			vArr3[i]=__vcorrcoef(vArr1,vArr2, mLen);
		}

		if(axis==0){
			free(vArr1);
			free(vArr2);
		}
	}
}

static void __mREZ(float *mArr1,int nLength,int mLength,int axis,float *vArr3,float (*func)(float *,int )){
	int nLen=0;
	int mLen=0;

	float *vArr=NULL;

	if(!vArr3){
		return;
	}

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;

			vArr=__vnew(mLen, NULL);
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			if(axis==0){
				for(int j=0;j<mLen;j++){
					if(axis==0){ // col遍历
						vArr[j]=mArr1[i+j*mLength];
					}
				}
			}
			else{
				vArr=mArr1+i*mLength;
			}

			vArr3[i]=func(vArr, mLen);
		}

		if(axis==0){
			free(vArr);
		}
	}
}

void __mrms(float *mArr1,int nLength,int mLength,int axis,float *vArr3){
	
	__mREZ(mArr1,nLength,mLength,axis,vArr3,__vrms);
}

void __menergy(float *mArr1,int nLength,int mLength,int axis,float *vArr3){

	__mREZ(mArr1,nLength,mLength,axis,vArr3,__venergy);
}

void __mzcr(float *mArr1,int nLength,int mLength,int axis,float *vArr3){

	__mREZ(mArr1,nLength,mLength,axis,vArr3,__vzcr);
}

void __munwrap(float *mArr1,int nLength,int mLength,int axis){
	int nLen=0;
	int mLen=0;

	float **vArr1=NULL;

	if(axis==0||axis==1){
		if(axis==0){
			nLen=mLength;
			mLen=nLength;

			vArr1=(float **)calloc(mLen, sizeof(float *));
		}
		else{
			nLen=nLength;
			mLen=mLength;
		}

		for(int i=0;i<nLen;i++){
			if(axis==0){
				for(int j=0;j<mLen;j++){
					if(axis==0){ // col遍历
						vArr1[j]=mArr1+(i+j*mLength);
					}
				}

				__vunwrap1(vArr1, mLen);
			}
			else{
				__vunwrap(mArr1+i*mLength, mLen, NULL);
			}
		}

		if(axis==0){
			free(vArr1);
		}
	}
}

float __vdot(float *vArr1,float *vArr2,int length){
	float value=0;

	for(int i=0;i<length;i++){
		value+=vArr1[i]*vArr2[i];
	}

	return value;
}

void __vmul(float *vArr1,float *vArr2,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]*vArr2[i];
	}
}

void __vdiv(float *vArr1,float *vArr2,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]/vArr2[i];
	}
}

void __vadd(float *vArr1,float *vArr2,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]+vArr2[i];
	}
}

void __vsub(float *vArr1,float *vArr2,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]-vArr2[i];
	}
}

void __vmul_value(float *vArr1,float value,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]*value;
	}
}

void __vdiv_value(float *vArr1,float value,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]/value;
	}
}

void __vadd_value(float *vArr1,float value,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]+value;
	}
}

void __vsub_value(float *vArr1,float value,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]-value;
	}
}

float __vnorm(float *vArr1,int length){
	float value=0;

	for(int i=0;i<length;i++){
		value+=vArr1[i]*vArr1[i];
	}
	value=sqrtf(value);

	return value;
}

// type 0 p 1/2/3... 1 Inf 2 -Inf
void __mnormalize(float *mArr1,int nLength,int mLength,int axis,int type,float p,float *mArr3){
	float *arr=NULL;

	int nLen=0;
	int mLen=0;

	float *vArr=NULL;

	if(mArr3){
		arr=mArr3;
	}
	else{
		arr=mArr1;
	}

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	vArr=__vnew(nLen, NULL);

	// 1. cal vArr
	for(int i=0;i<nLen;i++){
		vArr[i]=0;
		for(int j=0;j<mLen;j++){
			float _value=0;

			if(axis==0){ // 1*m
				_value=mArr1[i+j*mLength];
			}
			else{ // n*1
				_value=mArr1[i*mLength+j];
			}

			_value=fabsf(_value);
			if(type==0){ // P-L范式
				if(p==1.0){
					vArr[i]+=_value;
				}
				else if(p==2.0){
					vArr[i]+=_value*_value;
				}
				else{
					vArr[i]+=powf(_value, p);
				}
			}
			else if(type==1){ // max 
				if(j==0){
					vArr[i]=_value;
				}
				else{
					if(vArr[i]<_value){
						vArr[i]=_value;
					}
				}
			}
			else{ // min
				if(j==0){
					vArr[i]=_value;
				}
				else{
					if(vArr[i]>_value){
						vArr[i]=_value;
					}
				}
			}
		}

		if(type==0){ // P-L
			if(p!=1.0){
				if(p==2.0){
					vArr[i]=sqrtf(vArr[i]);
				}
				else{
					vArr[i]=powf(vArr[i], 1/p);
				}
			}
		}
	}

	// 2. div
	for(int i=0;i<nLen;i++){
		if(vArr[i]!=0){
			for(int j=0;j<mLen;j++){
				if(axis==0){ // 1*m
					arr[i+j*mLength]=mArr1[i+j*mLength]/vArr[i];
				}
				else{ // n*1
					arr[i*mLength+j]=mArr1[i*mLength+j]/vArr[i];
				}
			}
		}
	}

	free(vArr);
}

void __mnorm(float *mArr1,int nLength,int mLength,int axis,int type,float p,float *vArr2){
	float *arr=NULL;

	int nLen=0;
	int mLen=0;

	float *vArr=NULL;

	arr=mArr1;
	vArr=vArr2;

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	// 1. cal vArr
	for(int i=0;i<nLen;i++){
		vArr[i]=0;
		for(int j=0;j<mLen;j++){
			float _value=0;

			if(axis==0){ // 1*m
				_value=mArr1[i+j*mLength];
			}
			else{ // n*1
				_value=mArr1[i*mLength+j];
			}

			_value=fabsf(_value);
			if(type==0){ // P-L范式
				if(p==1.0){
					vArr[i]+=_value;
				}
				else if(p==2.0){
					vArr[i]+=_value*_value;
				}
				else{
					vArr[i]+=powf(_value, p);
				}
			}
			else if(type==1){ // max 
				if(j==0){
					vArr[i]=_value;
				}
				else{
					if(vArr[i]<_value){
						vArr[i]=_value;
					}
				}
			}
			else{ // min
				if(j==0){
					vArr[i]=_value;
				}
				else{
					if(vArr[i]>_value){
						vArr[i]=_value;
					}
				}
			}
		}

		if(type==0){ // P-L
			if(p!=1.0){
				if(p==2.0){
					vArr[i]=sqrtf(vArr[i]);
				}
				else{
					vArr[i]=powf(vArr[i], 1/p);
				}
			}
		}
	}
}

// 向量归一化 type 0 p 1/2/3... 1 Inf 2 -Inf
void __vnormalize(float *vArr1,int length,int type,float p,float *vArr2){
	float *arr=NULL;

	float value=0;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		float _value=0;

		_value=vArr1[i];
		_value=fabsf(_value);
		if(type==0){ // P-L
			if(p==1.0){
				value+=_value;
			}
			else if(p==2.0){
				value+=_value*_value;
			}
			else{
				value+=powf(_value, p);
			}
		}
		else if(type==1){ // max
			if(i==0){
				value=_value;
			}
			else{
				if(value<_value){
					value=_value;
				}
			}
		}
		else{ // min
			if(i==0){
				value=_value;
			}
			else{
				if(value>_value){
					value=_value;
				}
			}
		}
	}

	if(type==0){ // P-L
		if(p!=1.0){
			if(p==2.0){
				value=sqrtf(value);
			}
			else{
				value=powf(value, 1/p);
			}
		}
	}

	if(value!=0){
		for(int i=0;i<length;i++){
			arr[i]=vArr1[i]/value;
		}
	}
}

// 各种scale
// 不适用稀疏数据
void __vminmaxscale(float *vArr1,int length,float *vArr2){
	float min=0;
	float max=0;

	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	__vmin(vArr1, length, &min);
	__vmax(vArr1, length, &max);
	if(max>min){
		for(int i=0;i<length;i++){
			arr[i]=(vArr1[i]-min)/(max-min);
		}
	}
}

// 通用数据标准化 均值0方差1 不适用用稀疏数据
void __vstandscale(float *vArr1,int length,int type,float *vArr2){
	float mean=0;
	float std=0;

	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	mean=__vmean(vArr1, length);
	std=__vstd(vArr1, length, type);
	if(std){
		for(int i=0;i<length;i++){
			arr[i]=(vArr1[i]-mean)/std;
		}
	}
}

// 使用稀疏矩阵
void __vmaxabsscale(float *vArr1,int length,float *vArr2){
	float max=0;

	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	__vmaxabs(vArr1, length, &max);
	if(max){
		for(int i=0;i<length;i++){
			arr[i]=vArr1[i]/max;
		}
	}
}

void __vrobustscale(float *vArr1,int length,float *vArr2){
	float q1=0;
	float q3=0;
	float q2=0;

	int index=0;
	int mod=0;

	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	index=(length+1)/2-1;
	mod=(length+1)%2;
	if(index>=0){
		if(!mod){ // 整除
			q2=vArr1[index];
		}
		else{
			q2=(vArr1[index]+vArr1[index+1])/2;
		}
	}

	index=(length+1)/4-1;
	mod=(length+1)%4;
	if(index>=0){
		if(!mod){ // 整除
			q1=vArr1[index];
		}
		else{
			q1=(vArr1[index]+vArr1[index+1])/2;
		}
	}

	index=(length+1)*3/4-1;
	mod=(length+1)*3%4;
	if(index>=0){
		if(!mod){ // 整除
			q3=vArr1[index];
		}
		else{
			q3=(vArr1[index]+vArr1[index+1])/2;
		}
	}

	if(q3>q1){
		for(int i=0;i<length;i++){
			arr[i]=(vArr1[i]-q2)/(q3-q1);
		}
	}
}

void __vcenterscale(float *vArr1,int length,float *vArr2){
	float mean=0;

	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	mean=__vmean(vArr1, length);
	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]-mean;
	}
}

void __vmeanscale(float *vArr1,int length,float *vArr2){
	float min=0;
	float max=0;

	float mean=0;

	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	__vmin(vArr1, length, &min);
	__vmax(vArr1, length, &max);

	mean=__vmean(vArr1, length);
	if(max>min){
		for(int i=0;i<length;i++){
			arr[i]=(vArr1[i]-mean)/(max-min);
		}
	}
}

void __varctanscale(float *vArr1,int length,float *vArr2){
	float *arr=NULL;

	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=atanf(vArr1[i])/(M_PI/2);
	}
}

float __vsum(float *vArr1,int length){
	double value=0;

	for(int i=0;i<length;i++){
		value+=vArr1[i];
	}
	
	return value;
}

int __vsumi(int *vArr1,int length){
	int value=0;

	for(int i=0;i<length;i++){
		value+=vArr1[i];
	}
	
	return value;
}

int __vmin(float *vArr1,int length,float *value){
	int index=0;
	float min=0;

	if(!vArr1||!length){
		return -1;
	}

	min=vArr1[0];
	for(int i=1;i<length;i++){
		if(min>vArr1[i]){
			min=vArr1[i];
			index=i;
		}
	}

	if(value){
		*value=min;
	}

	return index;
}

int __vmax(float *vArr1,int length,float *value){
	int index=0;
	float max=0;

	if(!vArr1||!length){
		return -1;
	}

	max=vArr1[0];
	for(int i=1;i<length;i++){
		if(max<vArr1[i]){
			max=vArr1[i];
			index=i;
		}
	}

	if(value){
		*value=max;
	}

	return index;
}

int __vmini(int *vArr1,int length,int *value){
	int index=0;
	int min=0;

	if(!vArr1||!length){
		return -1;
	}

	min=vArr1[0];
	for(int i=1;i<length;i++){
		if(min>vArr1[i]){
			min=vArr1[i];
			index=i;
		}
	}

	if(value){
		*value=min;
	}

	return index;
}

int __vmaxi(int *vArr1,int length,int *value){
	int index=0;
	int max=0;

	if(!vArr1||!length){
		return -1;
	}

	max=vArr1[0];
	for(int i=1;i<length;i++){
		if(max<vArr1[i]){
			max=vArr1[i];
			index=i;
		}
	}

	if(value){
		*value=max;
	}

	return index;
}

int __vminabs(float *vArr1,int length,float *value){
	int index=0;
	float min=0;

	if(!vArr1||!length){
		return -1;
	}

	min=fabs(vArr1[0]);
	for(int i=1;i<length;i++){
		if(min>fabs(vArr1[i])){
			min=fabs(vArr1[i]);
			index=i;
		}
	}

	if(value){
		*value=min;
	}

	return index;
}

int __vmaxabs(float *vArr1,int length,float *value){
	int index=0;
	float max=0;

	if(!vArr1||!length){
		return -1;
	}

	max=fabs(vArr1[0]);
	for(int i=1;i<length;i++){
		if(max<fabs(vArr1[i])){
			max=fabs(vArr1[i]);
			index=i;
		}
	}

	if(value){
		*value=max;
	}

	return index;
}

float __vmean(float *vArr1,int length){
	float value=0;

	for(int i=0;i<length;i++){
		value+=vArr1[i];
	}
	value/=length;

	return value;
}

float __vmedian(float *vArr1,int length){
	float value=0;
	float *vArr2=NULL;

	int start=0;
	int len=0;

	if(!vArr1||!length){
		return 0;
	}

	if(length==1){
		return vArr1[0];
	}
	
	if(length==2){
		return (vArr1[0]+vArr1[1])/2;
	}

	vArr2=__vnew(length, NULL);	

	__vsort(vArr1, length, 0, vArr2);

	if(length&1){ // odd
		start=length/2;
		len=1;
	}
	else{
		start=length/2-1;
		len=2;
	}

	value=__vmean(vArr2+start, len);

	free(vArr2);
	return value;
}

float __vvar(float *vArr1,int length,int type){
	float mValue=0;
	float value=0;

	mValue=__vmean(vArr1, length);
	for(int i=0;i<length;i++){
		value+=(vArr1[i]-mValue)*(vArr1[i]-mValue);
	}

	value=value/(!type?length-1:length);
	return value;
}

float __vstd(float *vArr1,int length,int type){
	float value=0;

	value=__vvar(vArr1,length,type);
	value=sqrtf(value);
	return value;
}

float __vcov(float *vArr1,float *vArr2,int length,int type){
	float m1=0;
	float m2=0;

	float value=0;

	m1=__vmean(vArr1, length);
	m2=__vmean(vArr2, length);
	for(int i=0;i<length;i++){
		value+=(vArr1[i]-m1)*(vArr2[i]-m2);
	}

	value=value/(!type?length-1:length);
	return value;
}

float __vcorrcoef(float *vArr1,float *vArr2,int length){
	float c1=0;
	float s1=0;
	float s2=0;

	float value=0;

	c1=__vcov(vArr1, vArr2, length, 0);
	s1=__vstd(vArr1, length, 0);
	s2=__vstd(vArr2, length, 0);

	value=c1/(s1*s2);
	return value;
}

float __vrms(float *vArr1,int length){
	float value=0;

	for(int i=0;i<length;i++){
		value+=vArr1[i]*vArr1[i];
	}

	value/=length;
	value=sqrtf(value);

	return value;
}

float __venergy(float *vArr1,int length){
	float value=0;

	for(int i=0;i<length;i++){
		value+=vArr1[i]*vArr1[i];
	}

	return value;
}

float __vzcr(float *vArr1,int length){
	float rate=0;
	int num=0;

	if(vArr1&&length){
		for(int i=1;i<length;i++){
			if(vArr1[i]*vArr1[i-1]<0){
				num++;
			}
		}

		rate=1.0*num/length;
	}
	
	return rate;
}

void __vunwrap(float *vArr1,int length,float *vArr2){
	float sub=0;

	float mod=0;
	int t=0;

	float *arr=NULL;

	if(vArr1&&length){
		if(vArr2){
			arr=vArr2;
		}
		else{
			arr=vArr1;
		}

		arr[0]=vArr1[0];
		for(int i=1;i<length;i++){
			sub=fabsf(vArr1[i]-arr[i-1]);
			if(sub<M_PI){
				arr[i]=vArr1[i];
			}
			else{
				t=floorf(sub/(2*M_PI));
				mod=sub-t*2*M_PI;
				if(mod>M_PI){
					t++;
				}

				if(vArr1[i]>vArr1[i-1]){
					arr[i]=vArr1[i]-t*2*M_PI;
				}
				else{
					arr[i]=vArr1[i]+t*2*M_PI;
				}
			}
		}
	}
}

static void __vunwrap1(float **vArr1,int length){
	float sub=0;

	float mod=0;
	int t=0;

	if(vArr1&&length){
		for(int i=1;i<length;i++){
			sub=fabsf(vArr1[i][0]-vArr1[i-1][0]);
			if(sub>=M_PI){
				t=floorf(sub/(2*M_PI));
				mod=sub-t*2*M_PI;
				if(mod>M_PI){
					t++;
				}

				if(vArr1[i][0]>vArr1[i-1][0]){
					vArr1[i][0]=vArr1[i][0]-t*2*M_PI;
				}
				else{
					vArr1[i][0]=vArr1[i][0]+t*2*M_PI;
				}
			}
		}
	}
}

// type 0 升 1降
void __vsort(float *vArr1,int length,int type,float *vArr2){
	float *arr=NULL;

	if(vArr2&&vArr2!=vArr1){
		arr=vArr2;
		memcpy(arr, vArr1, sizeof(float )*length);
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		for(int j=i+1;j<length;j++){
			float _value=0;

			if(!type){ // 升
				if(arr[i]>arr[j]){
					_value=arr[i];
					arr[i]=arr[j];
					arr[j]=_value;
				}
			}
			else{ // 降
				if(arr[i]<arr[j]){
					_value=arr[i];
					arr[i]=arr[j];
					arr[j]=_value;
				}
			}
		}
	}
}

void __vsorti(int *vArr1,int length,int type,int *vArr2){
	int *arr=NULL;

	if(vArr2&&vArr2!=vArr1){
		arr=vArr2;
		memcpy(arr, vArr1, sizeof(int )*length);
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		for(int j=i+1;j<length;j++){
			int _value=0;

			if(!type){ // 升
				if(arr[i]>arr[j]){
					_value=arr[i];
					arr[i]=arr[j];
					arr[j]=_value;
				}
			}
			else{ // 降
				if(arr[i]<arr[j]){
					_value=arr[i];
					arr[i]=arr[j];
					arr[j]=_value;
				}
			}
		}
	}
}

// element math相关
void __vmap(float *vArr1,int length,void *callback,float *vArr2){
	float *arr=NULL;
	UniFunc call=NULL;

	call=(UniFunc )callback;
	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=call(vArr1[i]);
	}
}

void __vmap1(float *vArr1,int length,void *callback1,float value,float *vArr2){
	float *arr=NULL;
	UniFunc1 call=NULL;

	call=(UniFunc1 )callback1;
	if(vArr2){
		arr=vArr2;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=call(vArr1[i],value);
	}
}

// index/range/liespace/transpose/new相关
float *__vnew(int length,float *value){
	float *arr=NULL;
	float _value=0;

	arr=(float *)calloc(length, sizeof(float ));
	if(value){
		_value=*value;
	}

	if(_value!=0){
		for(int i=0;i<length;i++){
			arr[i]=_value;
		}
	}
	
	return arr;
}

// type 0 has stop 1 no stop
float *__vlinspace(float start,float stop,int length,int type){
	float *arr=NULL;
	float step=0;

	arr=__vnew(length,NULL);
	if(type==0){
		step=(stop-start)/(length-1>0?length-1:1);
	}
	else{
		step=(stop-start)/length;
	}
	
	for(int i=0;i<length;i++){
		arr[i]=start+i*step;
	}

	return arr;
}

int __varange(float start,float stop,float step,float **outArr){
	float *arr=NULL;
	int length=0;

	if(stop<=start||!outArr){
		return 0;
	}

	length=ceilf((stop-start)/step);
	arr=__vnew(length, NULL);
	for(int i=0;i<length;i++){
		arr[i]=start+i*step;
	}

	*outArr=arr;
	return length;
}

void __vfill(float *arr,int length,float value){

	for(int i=0;i<length;i++){
		arr[i]=value;
	}
}

void __vcopy(float *dstArr,float *srcArr,int length){

	memcpy(dstArr, srcArr, sizeof(float )*length);
}

void __mtrans(float *mArr1,int nLength,int mLength,float *mArr3){
	int index=0;

	if(!mArr3||mArr3==mArr1){
		return;
	}

	for(int i=0;i<mLength;i++){
		for(int j=0;j<nLength;j++){
			mArr3[index]=mArr1[i+j*mLength];
			index++;
		}
	}
}

float *__vget(float *vArr1,int length,VSlice *slice){
	float *vArr2=NULL;

	int nLen=0;

	int *rIndexArr=NULL;

	int rIndex=0;
	int index=0;

	if(!slice){
		return NULL;
	}
	
	if(slice->nLength>0){
		rIndexArr=slice->indexArr;
		nLen=slice->nLength;
	}
	else{
		nLen=length;
	}

	vArr2=__vnew(nLen, NULL);
	for(int i=0;i<nLen;i++){
		if(slice->nLength>0){
			rIndex=rIndexArr[i];
		}
		else{
			if(slice->nLength==0){  
				rIndex=i;
			}
			else{ // -1
				rIndex=nLen-1-i;
			}
		}

		vArr2[index]=vArr1[rIndex];
		index++;
	}

	return vArr2;
}

void __vset(float *vArr1,int length,VSlice *slice,float *vArr2){
	int nLen=0;

	int *rIndexArr=NULL;

	int rIndex=0;
	int index=0;

	if(!slice){
		return ;
	}
	
	if(slice->nLength>0){
		rIndexArr=slice->indexArr;
		nLen=slice->nLength;
	}
	else{
		nLen=length;
	}

	for(int i=0;i<nLen;i++){
		if(slice->nLength>0){
			rIndex=rIndexArr[i];
		}
		else{
			if(slice->nLength==0){  
				rIndex=i;
			}
			else{ // -1
				rIndex=nLen-1-i;
			}
		}

		vArr1[rIndex]=vArr2[index];
		index++;
	}
}

/***
	[:,:]
	nLength => i row
	mLength => j col
	0 => range(0,nLength/mLength)
	-1 => range(nLength/mLength,0)
****/
float *__mget(float *mArr1,int nLength,int mLength,VSlice *slice){
	float *mArr2=NULL;

	int nLen=0;
	int mLen=0;

	int *rIndexArr=NULL;
	int *cIndexArr=NULL;

	int rIndex=0;
	int cIndex=0;

	int index=0;

	if(!slice){
		return NULL;
	}
	
	if(slice->nLength>0){
		rIndexArr=slice->indexArr;
		nLen=slice->nLength;
	}
	else{
		nLen=nLength;
	}

	if(slice->mLength>0){
		if(slice->nLength>0){
			cIndexArr=slice->indexArr+slice->nLength;
		}
		else{
			cIndexArr=slice->indexArr;
		}
		mLen=slice->mLength;
	}
	else{
		mLen=mLength;
	}

	mArr2=__vnew(nLen*mLen, NULL);
	for(int i=0;i<nLen;i++){
		if(slice->nLength>0){
			rIndex=rIndexArr[i];
		}
		else{
			if(slice->nLength==0){  
				rIndex=i;
			}
			else{ // -1
				rIndex=nLen-1-i;
			}
		}

		for(int j=0;j<mLen;j++){
			if(slice->mLength>0){
				cIndex=cIndexArr[j];
			}
			else{
				if(slice->mLength==0){  
					cIndex=j;
				}
				else{ // -1
					cIndex=mLen-1-j;
				}
			}

			mArr2[index]=mArr1[rIndex*mLength+cIndex];
			index++;
		}
	}

	return mArr2;
}

void __mset(float *mArr1,int nLength,int mLength,VSlice *slice,float *mArr2){
	int nLen=0;
	int mLen=0;

	int *rIndexArr=NULL;
	int *cIndexArr=NULL;

	int rIndex=0;
	int cIndex=0;

	int index=0;

	if(!slice){
		return ;
	}
	
	if(slice->nLength>0){
		rIndexArr=slice->indexArr;
		nLen=slice->nLength;
	}
	else{
		nLen=nLength;
	}

	if(slice->mLength>0){
		if(slice->nLength>0){
			cIndexArr=slice->indexArr+slice->nLength;
		}
		else{
			cIndexArr=slice->indexArr;
		}
		mLen=slice->mLength;
	}
	else{
		mLen=mLength;
	}

	for(int i=0;i<nLen;i++){
		if(slice->nLength>0){
			rIndex=rIndexArr[i];
		}
		else{
			if(slice->nLength==0){  
				rIndex=i;
			}
			else{ // -1
				rIndex=nLen-1-i;
			}
		}

		for(int j=0;j<mLen;j++){
			if(slice->mLength>0){
				cIndex=cIndexArr[j];
			}
			else{
				if(slice->mLength==0){  
					cIndex=j;
				}
				else{ // -1
					cIndex=mLen-1-j;
				}
			}

			mArr1[rIndex*mLength+cIndex]=mArr2[index];
			index++;
		}
	}
}

void __mset_value(float *mArr1,int nLength,int mLength,VSlice *slice,float value){
	int nLen=0;
	int mLen=0;

	int *rIndexArr=NULL;
	int *cIndexArr=NULL;

	int rIndex=0;
	int cIndex=0;

	if(!slice){
		return ;
	}
	
	if(slice->nLength>0){
		rIndexArr=slice->indexArr;
		nLen=slice->nLength;
	}
	else{
		nLen=nLength;
	}

	if(slice->mLength>0){
		if(slice->nLength>0){
			cIndexArr=slice->indexArr+slice->nLength;
		}
		else{
			cIndexArr=slice->indexArr;
		}
		mLen=slice->mLength;
	}
	else{
		mLen=mLength;
	}

	for(int i=0;i<nLen;i++){
		if(slice->nLength>0){
			rIndex=rIndexArr[i];
		}
		else{
			if(slice->nLength==0){  
				rIndex=i;
			}
			else{ // -1
				rIndex=nLen-1-i;
			}
		}

		for(int j=0;j<mLen;j++){
			if(slice->mLength>0){
				cIndex=cIndexArr[j];
			}
			else{
				if(slice->mLength==0){  
					cIndex=j;
				}
				else{ // -1
					cIndex=mLen-1-j;
				}
			}

			mArr1[rIndex*mLength+cIndex]=value;
		}
	}
}

void __mset_vector(float *mArr1,int nLength,int mLength,VSlice *slice,int axis,float *vArr1){
	int nLen=0;
	int mLen=0;

	int nLength1=0;
	int mLength1=0;

	int *rIndexArr=NULL;
	int *cIndexArr=NULL;

	int rIndex=0;
	int cIndex=0;

	if(!slice){
		return ;
	}
	
	if(slice->nLength>0){
		rIndexArr=slice->indexArr;
		nLength1=slice->nLength;
	}
	else{
		nLength1=nLength;
	}

	if(slice->mLength>0){
		if(slice->nLength>0){
			cIndexArr=slice->indexArr+slice->nLength;
		}
		else{
			cIndexArr=slice->indexArr;
		}
		mLength1=slice->mLength;
	}
	else{
		mLength1=mLength;
	}

	if(axis==0){
		nLen=mLength1;
		mLen=nLength1;
	}
	else{
		nLen=nLength1;
		mLen=mLength1;
	}

	for(int i=0;i<nLen;i++){
		if(axis==0){ // ==mLength1
			if(slice->mLength>0){
				rIndex=cIndexArr[i];
			}
			else{
				if(slice->mLength==0){  
					rIndex=i;
				}
				else{ // -1
					rIndex=mLength1-1-i;
				}
			}
		}
		else{ // ==nLength1
			if(slice->nLength>0){
				rIndex=rIndexArr[i];
			}
			else{
				if(slice->nLength==0){  
					rIndex=i;
				}
				else{ // -1
					rIndex=nLength1-1-i;
				}
			}
		}

		for(int j=0;j<mLen;j++){
			if(axis==0){ // ==nLength1
				if(slice->nLength>0){
					cIndex=rIndexArr[j];
				}
				else{
					if(slice->nLength==0){  
						cIndex=j;
					}
					else{ // -1
						cIndex=nLength1-1-j;
					}
				}
			}
			else{ // ==mLength1
				if(slice->mLength>0){
					cIndex=cIndexArr[j];
				}
				else{
					if(slice->mLength==0){  
						cIndex=j;
					}
					else{ // -1
						cIndex=mLength1-1-j;
					}
				}
			}

			if(axis==0){
				mArr1[rIndex+cIndex*mLength]=vArr1[i];
			}
			else{
				mArr1[rIndex*mLength+cIndex]=vArr1[i];
			}
		}
	}
}

// astype
void __vf2i(float *vArr1,int length,int *vArr2){
	
	if(!vArr1||!vArr2){
		return;
	}

	for(int i=0;i<length;i++){
		vArr2[i]=floorf(vArr1[i]);
	}
}

void __vi2f(int *vArr1,int length,float *vArr2){
	
	if(!vArr1||!vArr2){
		return;
	}

	for(int i=0;i<length;i++){
		vArr2[i]=vArr1[i];
	}
}

// repeat axis 0 1*m 1 n*1
void __vrepeat(float *vArr1,int length,int num,float *vArr3){

	// if(axis==0){ // 1*m => num*length
	// 	for(int i=0;i<num;i++){
	// 		for(int j=0;j<length;j++){
	// 			mArr3[i*length+j]=vArr1[j];
	// 		}
	// 	}
	// }
	// else{ // n*1 => length*num
	// 	for(int i=0;i<num;i++){
	// 		for(int j=0;j<length;j++){
	// 			mArr3[i+j*num]=vArr1[j];
	// 		}
	// 	}
	// }

	for(int i=0;i<num;i++){
		for(int j=0;j<length;j++){
			vArr3[i*length+j]=vArr1[j];
		}
	}
}

void __mrepeat(float *mArr1,int nLength,int mLength,int axis,int num,float *mArr3){
	int length=0;

	if(axis==0){ // n*m => (n*num)*m
		length=nLength*num;
		for(int i=0;i<length;i++){
			for(int j=0;j<mLength;j++){
				mArr3[i*mLength+j]=mArr1[i*mLength%(nLength*mLength)+j];
			}
		}
	}
	else{ // n*m => n*(m*num)
		length=mLength*num;
		for(int i=0;i<length;i++){
			for(int j=0;j<nLength;j++){
				mArr3[i+j*length]=mArr1[i%mLength+j*mLength];
			}
		}
	}
}

void __vconcat(float *vArr1,float *vArr2,int length1,int length2,float *vArr3){

	for(int i=0;i<length1;i++){
		vArr3[i]=vArr1[i];
	}

	for(int i=0;i<length2;i++){
		vArr3[i+length1]=vArr2[i];
	}
}

void __mconcat(float *mArr1,float *mArr2,
				int nLength1,int mLength1,
				int nLength2,int mLength2,
				int axis,
				float *mArr3){
	if(axis==0){ // n1*m1 n2*m2 m1==m2
		if(mLength1!=mLength2){
			return;
		}

		for(int i=0;i<nLength1;i++){
			for(int j=0;j<mLength1;j++){
				mArr3[i*mLength1+j]=mArr1[i*mLength1+j];
			}
		}

		for(int i=0;i<nLength2;i++){
			for(int j=0;j<mLength2;j++){
				mArr3[(i+nLength1)*mLength2+j]=mArr2[i*mLength2+j];
			}
		}
	}
	else{ // n1*m1 n2*m2 n1==n2
		if(nLength1!=nLength2){
			return;
		}
		
		for(int i=0;i<mLength1;i++){
			for(int j=0;j<nLength1;j++){
				mArr3[i+j*(mLength1+mLength2)]=mArr1[i+j*mLength1];
			}
		}

		for(int i=0;i<mLength2;i++){
			for(int j=0;j<nLength2;j++){
				mArr3[i+mLength1+j*(mLength1+mLength2)]=mArr2[i+j*mLength2];
			}
		}	
	}
}

void __mcut(float *mArr1,int nLength,int mLength,
			int nIndex,int nLength1,
			int mIndex,int mLength1,
			float *mArr3){
	float *arr=NULL;
	int index=0;

	if(nIndex<=0||
		nIndex>nLength-1||
		mIndex<=0||
		mIndex>mLength-1||
		nIndex+nLength1>nLength||
		mIndex+mLength1>mLength){
		return;
	}

	if(mArr3){
		arr=mArr3;
	}
	else{
		arr=mArr1;
	}

	for(int i=nIndex;i<nLength1+nIndex;i++){
		for(int j=mIndex;j<mLength1+mIndex;j++){
			arr[index]=mArr1[i*mLength+j];
			index++;
		}
	}
}

// debug vector/matrix相关
void __vdebug(float *vArr1,int length,int type){
	if(type){
		printf("vector is:\n");

		printf("	");
		for(int i=0;i<length;i++){
			printf("%d:%f ,",i,vArr1[i]);
		}
	}
	else{
		printf("vector([");
		for(int i=0;i<length-1;i++){
			printf("%f, ",vArr1[i]);
		}
		printf("%f",vArr1[length-1]);
		printf("])\n");
	}
}

void __mdebug(float *mArr1,int nLength,int mLength,int type){
	if(type){
		printf("matrix is:\n");

		for(int i=0;i<nLength;i++){
			printf("	%d:\n",i);
			printf("		");
			for(int j=0;j<mLength;j++){
				printf("%d:%f ,",j,mArr1[i*mLength+j]);
			}

			printf("\n\n");
		}
	}
	else{
		printf("matrix([");
		for(int i=0;i<nLength-1;i++){
			printf("[");

			for(int j=0;j<mLength-1;j++){
				printf("%f ,",mArr1[i*mLength+j]);
			}
			printf("%f",mArr1[i*mLength+mLength-1]);

			printf("],");
			printf("\n        ");
		}

		printf("[");
		for(int j=0;j<mLength-1;j++){
			printf("%f ,",mArr1[(nLength-1)*mLength+j]);
		}
		printf("%f",mArr1[(nLength-1)*mLength+mLength-1]);

		printf("]");
		printf("])\n");
	}
}

void __mdiff(float *mArr1,int nLength,int mLength,int axis,int *order,float *mArr3){
	int nLen=0;
	int mLen=0;

	int _order=1;

	if(!mArr3){
		return;
	}

	if(order){
		_order=*order;
	}

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	for(int i=0;i<nLength*mLength;i++){
		mArr3[i]=mArr1[i];
	}

	while(_order>0&&mLen>1){
		for(int i=0;i<nLen;i++){
			for(int j=1;j<mLen;j++){
				if(axis==0){ // 1*m
					mArr3[i+(j-1)*mLength]=mArr3[i+j*mLength]-mArr3[i+(j-1)*mLength];
				}
				else{ // n*1
					mArr3[i*mLength+j-1]=mArr3[i*mLength+j]-mArr3[i*mLength+j-1];
					mLength--;
				}
			}
		}

		_order--;
		mLen--;
	}
}

// step 1
void __mdiff2(float *mArr1,int nLength,int mLength,int axis,int *step,float *mArr3){
	int nLen=0;
	int mLen=0;

	int _step=1;

	if(!mArr3){
		return;
	}

	if(step){
		_step=*step;
	}

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	for(int i=0;i<nLength*mLength;i++){
		mArr3[i]=mArr1[i];
	}

	if(_step>0){
		for(int i=0;i<nLen;i++){
			for(int j=0;j<mLen;j++){
				if(axis==0){ // 1*m
					if(j<_step){
						mArr3[i+j*mLength]=0;
					}
					else{
						mArr3[i+j*mLength]=mArr3[i+j*mLength]-mArr1[i+(j-_step)*mLength];
					}
				}
				else{ // n*1
					if(j<_step){
						mArr3[i*mLength+j]=0;
					}
					else{
						mArr3[i*mLength+j]=mArr3[i*mLength+j]-mArr1[i*mLength+j-_step];
					}
				}
			}
		}
	}
}

void __mdiff3(float *mArr1,float *mArr2,int nLength,int mLength,int axis,int *step,float *mArr3){


}

// type 0 max 1 median
static void __mxfilter(float *mArr1,int nLength,int mLength,int type,int axis,int order,float *mArr3){
	int nLen=0;
	int mLen=0;

	float *vArr1=NULL;
	float *vArr2=NULL;

	if(!mArr3||order<1){
		return;
	}

	if(axis==0){
		nLen=mLength;
		mLen=nLength;
	}
	else{
		nLen=nLength;
		mLen=mLength;
	}

	vArr1=__vnew(mLen, NULL);
	vArr2=__vnew(mLen, NULL);

	for(int i=0;i<nLen;i++){
		if(axis==0){ // 1*m
			for(int j=0;j<mLen;j++){
				vArr1[j]=mArr1[i+j*mLength];
			}

			if(type==0){
				__vmaxfilter(vArr1,mLen,order,vArr2);
			}
			else{
				__vmedianfilter(vArr1,mLen,order,vArr2);
			}

			for(int j=0;j<mLen;j++){
				mArr3[i+j*mLength]=vArr2[j];
			}
		}
		else{
			if(type==0){
				__vmaxfilter(mArr1+i*mLength,mLen,order,vArr2);
			}
			else{
				__vmedianfilter(mArr1+i*mLength,mLen,order,vArr2);
			}

			for(int j=0;j<mLen;j++){
				mArr3[i*mLength+j]=vArr2[j];
			}
		}
	}

	free(vArr1);
	free(vArr2);
}

void __mmaxfilter(float *mArr1,int nLength,int mLength,int axis,int order,float *mArr3){
	
	__mxfilter(mArr1,nLength,mLength,0,axis,order,mArr3);
}

void __mmedianfilter(float *mArr1,int nLength,int mLength,int axis,int order,float *mArr3){
	
	if(order>1&&(order&1)){
		__mxfilter(mArr1,nLength,mLength,1,axis,order,mArr3);
	}
}

void __vdiff(float *vArr1,int length,int *order,float *vArr3){
	int _order=1;

	if(!vArr3){
		return;
	}

	if(order){
		_order=*order;
	}

	for(int i=0;i<length;i++){
		vArr3[i]=vArr1[i];
	}

	while(_order>0&&length>1){
		for(int i=1;i<length;i++){
			vArr3[i-1]=vArr3[i]-vArr3[i-1];
		}

		_order--;
		length--;
	}
}

// step 1
void __vdiff2(float *vArr1,int length,int *step,float *vArr3){
	int _step=1;

	if(!vArr3){
		return;
	}

	if(step){
		_step=*step;
	}

	if(_step>0){
		memset(vArr3, 0, sizeof(float )*_step);
		for(int i=_step;i<length;i++){
			vArr3[i]=vArr1[i];
		}

		for(int i=_step;i<length;i++){
			vArr3[i]=vArr3[i]-vArr1[i-_step];
		}
	}
}

void __vmaxfilter(float *vArr1,int length,int order,float *vArr3){
	int left=0;
	int right=0;

	int start=0;
	int end=0;

	if(order<1||!vArr3){
		return;
	}

	left=order/2;
	right=order-left;
	for(int i=0;i<length;i++){
		start=(i-left>=0?i-left:0);
		end=(i-1+right<=length-1?i-1+right:length-1);
		vArr3[i]=__arr_max(vArr1+start, end-start+1);
	}
}

void __vmedianfilter(float *vArr1,int length,int order,float *vArr3){
	float *arr1=NULL;

	int len=0;

	if((order&1)==0||order<2||!vArr3){
		return;
	}

	len=order/2;
	arr1=__vnew(length+2*len, NULL);
	memcpy(arr1+len, vArr1, sizeof(float )*length);

	for(int i=len,j=0;i<length+len;i++,j++){
		float value=0;

		value=__vmedian(arr1+(i-len), order);
		vArr3[j]=value;
	}

	free(arr1);
}

static float __arr_max(float *vArr,int length){
	float max=0;

	if(!vArr||!length){
		return 0;
	}

	max=vArr[0];
	for(int i=1;i<length;i++){
		if(max<vArr[i]){
			max=vArr[i];
		}
	}

	return max;
}




