// 

#include <string.h>
#include <math.h>

#include "flux_complex.h"

// 针对 type 0即符合乘法
void __mcdot(float *mRealArr1,float *mImageArr1,
			float *mRealArr2,float *mImageArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			float *mRealArr3,float *mImageArr3){
	if(mLength1!=nLength2){
		return;
	}

	for(int i=0;i<nLength1;i++){
		for(int j=0;j<mLength2;j++){
			double _value1=0;
			double _value2=0;

			float r1=0;
			float r2=0;

			float i1=0;
			float i2=0;

			for(int k=0;k<mLength1;k++){ // arr1.row*arr2.col
				r1=mRealArr1[i*mLength1+k];
				r2=mRealArr2[j+k*mLength2];

				i1=mImageArr1[i*mLength1+k];
				i2=mImageArr2[j+k*mLength2];

				_value1+=(r1*r2-i1*i2);
				_value2+=(i1*r2+r1*i2);
			}

			mRealArr3[i*mLength2+j]=_value1;
			mImageArr3[i*mLength2+j]=_value2;
		}
	}
}

// 针对 type 1即m2转置符合乘法
void __mcdot1(float *mRealArr1,float *mImageArr1,
			float *mRealArr2,float *mImageArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			float *mRealArr3,float *mImageArr3){
	if(mLength1!=mLength2){
		return;
	}

	for(int i=0;i<nLength1;i++){
		for(int j=0;j<nLength2;j++){
			double _value1=0;
			double _value2=0;

			float r1=0;
			float r2=0;

			float i1=0;
			float i2=0;

			for(int k=0;k<mLength1;k++){ // arr1.row*arr2.row
				r1=mRealArr1[i*mLength1+k];
				r2=mRealArr2[j*mLength2+k];

				i1=mImageArr1[i*mLength1+k];
				i2=mImageArr2[j*mLength2+k];

				_value1+=(r1*r2-i1*i2);
				_value2+=(i1*r2+r1*i2);
			}

			mRealArr3[i*nLength2+j]=_value1;
			mImageArr3[i*nLength2+j]=_value2;
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
int __mcdot2(float *mRealArr1,float *mImageArr1,
			float *mRealArr2,float *mImageArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			int *type,
			float *mRealArr3,float *mImageArr3){
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
			double _value1=0;
			double _value2=0;

			float r1=0;
			float r2=0;

			float i1=0;
			float i2=0;

			if(_type==0){ // n1*m2
				for(int k=0;k<mLength1;k++){ // arr1.row*arr2.col
					r1=mRealArr1[i*mLength1+k];
					r2=mRealArr2[j+k*mLength2];

					i1=mImageArr1[i*mLength1+k];
					i2=mImageArr2[j+k*mLength2];

					_value1+=(r1*r2-i1*i2);
					_value2+=(i1*r2+r1*i2);
				}
			}
			else if(_type==1){ // n1*n2 
				for(int k=0;k<mLength1;k++){ // arr1.row*arr2.T.col(arr2.row)
					r1=mRealArr1[i*mLength1+k];
					r2=mRealArr2[j*mLength2+k];

					i1=mImageArr1[i*mLength1+k];
					i2=mImageArr2[j*mLength2+k];

					_value1+=(r1*r2-i1*i2);
					_value2+=(i1*r2+r1*i2);
				}
			}
			else if(_type==2){ // m1*m2 
				for(int k=0;k<nLength1;k++){ // arr1.T.row(arr1.col)*arr2.col
					r1=mRealArr1[i+k*mLength1];
					r2=mRealArr2[j+k*mLength2];

					i1=mImageArr1[i+k*mLength1];
					i2=mImageArr2[j+k*mLength2];

					_value1+=(r1*r2-i1*i2);
					_value2+=(i1*r2+r1*i2);
				}
			}
			else{ // m1*n2
				for(int k=0;k<nLength1;k++){ // arr1.T.row(arr1.col)*arr2.T.col(arr2.row)
					r1=mRealArr1[i+k*mLength1];
					r2=mRealArr2[j*mLength2+k];

					i1=mImageArr1[i+k*mLength1];
					i2=mImageArr2[j*mLength2+k];

					_value1+=(r1*r2-i1*i2);
					_value2+=(i1*r2+r1*i2);
				}
			}

			mRealArr3[i*mLen3+j]=_value1;
			mImageArr3[i*mLen3+j]=_value2;
		}
	}

	return status;
}

void __mcdiv(float *mRealArr1,float *mImageArr1,
			float *mRealArr2,float *mImageArr2,
			int nLength,int mLength,
			float *mRealArr3,float *mImageArr3){
	float *rArr=NULL;
	float *iArr=NULL;

	float real1=0;
	float image1=0;

	float real2=0;
	float image2=0;

	float real3=0;
	float image3=0;

	if(mRealArr3&&mImageArr3){
		rArr=mRealArr3;
		iArr=mImageArr3;
	}
	else{
		rArr=mRealArr1;
		iArr=mImageArr1;
	}

	for(int i=0;i<nLength;i++){
		for(int j=0;j<mLength;j++){
			real1=mRealArr1[i*mLength+j];
			image1=mImageArr1[i*mLength+j];

			real2=mRealArr2[i*mLength+j];
			image2=mImageArr2[i*mLength+j];

			__complexDiv(real1,image1,real2,image2,&real3,&image3);
			rArr[i*mLength+j]=real3;
			iArr[i*mLength+j]=image3;
		}
	}
}

void __mccut(float *mRealArr1,float *mImageArr1,
			int nLength,int mLength,
			int nIndex,int nLength1,
			int mIndex,int mLength1,
			float *mRealArr3,float *mImageArr3){
	float *rArr=NULL;
	float *iArr=NULL;
	int index=0;

	if(nIndex<0||
		nIndex>nLength-1||
		mIndex<0||
		mIndex>mLength-1||
		nIndex+nLength1>nLength||
		mIndex+mLength1>mLength){
		return;
	}

	if(mRealArr3&&mImageArr3){
		rArr=mRealArr3;
		iArr=mImageArr3;
	}
	else{
		rArr=mRealArr1;
		iArr=mImageArr1;
	}

	for(int i=nIndex;i<nLength1+nIndex;i++){
		for(int j=mIndex;j<mLength1+mIndex;j++){
			rArr[index]=mRealArr1[i*mLength+j];
			iArr[index]=mImageArr1[i*mLength+j];
			index++;
		}
	}
}

void __vcnew(int length,float *value,float **realArr,float **imageArr){
	float *rArr=NULL;
	float *iArr=NULL;

	float _value=0;

	if(value){
		_value=*value;
	}

	if(realArr){
		rArr=(float *)calloc(length, sizeof(float ));
		for(int i=0;i<length;i++){
			rArr[i]=_value;
		}

		*realArr=rArr;
	}

	if(imageArr){
		iArr=(float *)calloc(length, sizeof(float ));
		for(int i=0;i<length;i++){
			iArr[i]=_value;
		}

		*imageArr=iArr;
	}
}

void __vcabs(float *realArr1,float *imageArr1,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=imageArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=sqrtf(realArr1[i]*realArr1[i]+imageArr1[i]*imageArr1[i]);
	}
}

void __vcsquare(float *realArr1,float *imageArr1,int length,float *vArr3){
	float *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=imageArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=realArr1[i]*realArr1[i]+imageArr1[i]*imageArr1[i];
	}
}

// n*m complex=>n*m1 
void __mcabs(float *mRealArr,float *mImageArr,int nLength,int mLength,int axis,float *mArr2){
	float *mArr=NULL;

	int nLen=0;
	int mLen=0;

	if(mArr2){
		mArr=mArr2;
	}
	else{
		mArr=mImageArr;
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
			if(axis==0){
				mArr2[i+j*mLen]=sqrtf(mRealArr[i+j*mLen]*mRealArr[i+j*mLen]+
										mImageArr[i+j*mLen]*mImageArr[i+j*mLen]);
			}
			else{
				mArr2[i*mLen+j]=sqrtf(mRealArr[i*mLen+j]*mRealArr[i*mLen+j]+
										mImageArr[i*mLen+j]*mImageArr[i*mLen+j]);
			}
		}
	}
}

void __mcabs1(float *mRealArr,float *mImageArr,int nLength,int mLength,int axis,int cutLength,float *mArr2){
	float *mArr=NULL;

	int nLen=0;
	int mLen=0;

	if(mArr2){
		mArr=mArr2;
	}
	else{
		mArr=mImageArr;
	}

	if(axis==0){
		if(cutLength<0||cutLength>nLength){
			return;
		}

		nLen=cutLength;
		mLen=nLength;
	}
	else{
		if(cutLength<0||cutLength>mLength){
			return;
		}

		nLen=nLength;
		mLen=cutLength;
	}

	for(int i=0;i<nLen;i++){
		for(int j=0;j<mLen;j++){
			if(axis==0){
				mArr2[i+j*mLen]=sqrtf(mRealArr[i+j*mLen]*mRealArr[i+j*mLen]+
										mImageArr[i+j*mLen]*mImageArr[i+j*mLen]);
			}
			else{
				mArr2[i*mLen+j]=sqrtf(mRealArr[i*mLength+j]*mRealArr[i*mLength+j]+
										mImageArr[i*mLength+j]*mImageArr[i*mLength+j]);
			}
		}
	}
}

// n*m complex=>n*m1 
void __mcabs2(float *mRealArr,float *mImageArr,int nLength,int mLength,int mLength2,float *mArr2){
	float *mArr=NULL;

	if(mArr2){
		mArr=mArr2;
	}
	else{
		mArr=mImageArr;
	}

	if(mLength2<0||mLength2>mLength){
		return;
	}

	for(int i=0;i<nLength;i++){
		for(int j=0;j<mLength2;j++){
			mArr2[i*mLength2+j]=sqrtf(mRealArr[i*mLength+j]*mRealArr[i*mLength+j]+
										mImageArr[i*mLength+j]*mImageArr[i*mLength+j]);
		}
	}

}

void __mcsquare(float *mRealArr,float *mImageArr,int nLength,int mLength,int axis,float *mArr2){
	float *mArr=NULL;

	int nLen=0;
	int mLen=0;

	if(mArr2){
		mArr=mArr2;
	}
	else{
		mArr=mImageArr;
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
			if(axis==0){
				mArr2[i+j*mLen]=mRealArr[i+j*mLen]*mRealArr[i+j*mLen]+
										mImageArr[i+j*mLen]*mImageArr[i+j*mLen];
			}
			else{
				mArr2[i*mLen+j]=mRealArr[i*mLen+j]*mRealArr[i*mLen+j]+
										mImageArr[i*mLen+j]*mImageArr[i*mLen+j];
			}
		}
	}
}

void __mcsquare1(float *mRealArr,float *mImageArr,int nLength,int mLength,int axis,int cutLength,float *mArr2){
	float *mArr=NULL;

	int nLen=0;
	int mLen=0;

	if(mArr2){
		mArr=mArr2;
	}
	else{
		mArr=mImageArr;
	}

	if(axis==0){
		if(cutLength<0||cutLength>nLength){
			return;
		}

		nLen=cutLength;
		mLen=nLength;
	}
	else{
		if(cutLength<0||cutLength>mLength){
			return;
		}

		nLen=nLength;
		mLen=cutLength;
	}

	for(int i=0;i<nLen;i++){
		for(int j=0;j<mLen;j++){
			if(axis==0){
				mArr2[i+j*mLen]=mRealArr[i+j*mLen]*mRealArr[i+j*mLen]+
										mImageArr[i+j*mLen]*mImageArr[i+j*mLen];
			}
			else{
				mArr2[i*mLen+j]=mRealArr[i*mLength+j]*mRealArr[i*mLength+j]+
										mImageArr[i*mLength+j]*mImageArr[i*mLength+j];
			}
		}
	}
}

void __mcsquare2(float *mRealArr,float *mImageArr,int nLength,int mLength,int mLength2,float *mArr2){
	float *mArr=NULL;

	if(mArr2){
		mArr=mArr2;
	}
	else{
		mArr=mImageArr;
	}

	if(mLength2<0||mLength2>mLength){
		return;
	}

	for(int i=0;i<nLength;i++){
		for(int j=0;j<mLength2;j++){
			mArr2[i*mLength2+j]=mRealArr[i*mLength+j]*mRealArr[i*mLength+j]+
										mImageArr[i*mLength+j]*mImageArr[i*mLength+j];
		}
	}
}

// add/sub/mul/div
void __vcmul(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *realArr3,float *imageArr3){

	float *rArr=NULL;
	float *iArr=NULL;

	if(!realArr3||!imageArr3){
		rArr=realArr1;
		iArr=imageArr1;
	}
	else{
		rArr=realArr3;
		iArr=imageArr3;
	}

	for(int i=0;i<length;i++){
		__complexMul(realArr1[i], imageArr1[i], realArr2[i], imageArr2[i], rArr+i, iArr+i);
	}
}

void __vcdiv(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *realArr3,float *imageArr3){
	float *rArr=NULL;
	float *iArr=NULL;

	if(!realArr3||!imageArr3){
		rArr=realArr1;
		iArr=imageArr1;
	}
	else{
		rArr=realArr3;
		iArr=imageArr3;
	}

	for(int i=0;i<length;i++){
		__complexDiv(realArr1[i], imageArr1[i], realArr2[i], imageArr2[i], rArr+i, iArr+i);
	}
}

void __vcadd(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *realArr3,float *imageArr3){
	float *rArr=NULL;
	float *iArr=NULL;

	if(!realArr3||!imageArr3){
		rArr=realArr1;
		iArr=imageArr1;
	}
	else{
		rArr=realArr3;
		iArr=imageArr3;
	}

	for(int i=0;i<length;i++){
		rArr[i]=realArr1[i]+realArr2[i];
		iArr[i]=imageArr1[i]+imageArr2[i];
	}
}

void __vcsub(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *realArr3,float *imageArr3){
	float *rArr=NULL;
	float *iArr=NULL;

	if(!realArr3||!imageArr3){
		rArr=realArr1;
		iArr=imageArr1;
	}
	else{
		rArr=realArr3;
		iArr=imageArr3;
	}

	for(int i=0;i<length;i++){
		rArr[i]=realArr1[i]-realArr2[i];
		iArr[i]=imageArr1[i]-imageArr2[i];
	}
}

// debug vector/matrix相关
void __vcdebug(float *realArr,float *imageArr,int length,int type){
	if(type){
		printf("vector complex is:\n");

		printf("	");
		for(int i=0;i<length;i++){
			printf("%d:%f+j%f ,",i,realArr[i],imageArr[i]);
		}
	}
	else{
		printf("vector([");
		for(int i=0;i<length-1;i++){
			printf("%f+j%f, ",realArr[i],imageArr[i]);
		}
		printf("%f+j%f",realArr[length-1],imageArr[length-1]);
		printf("])\n");
	}

}

void __mcdebug(float *realArr,float *imageArr,int nLength,int mLength,int type){
	if(type){
		printf("matrix is:\n");

		for(int i=0;i<nLength;i++){
			printf("	%d:\n",i);
			printf("		");
			for(int j=0;j<mLength;j++){
				printf("%d:%f+j%f ,",j,realArr[i*mLength+j],imageArr[i*mLength+j]);
			}

			printf("\n\n");
		}
	}
	else{
		printf("matrix([");
		for(int i=0;i<nLength-1;i++){
			printf("[");

			for(int j=0;j<mLength-1;j++){
				printf("%f+j%f ,",realArr[i*mLength+j],imageArr[i*mLength+j]);
			}
			printf("%f++j%f",realArr[i*mLength+mLength-1],imageArr[i*mLength+mLength-1]);

			printf("],");
			printf("\n        ");
		}

		printf("[");
		for(int j=0;j<mLength-1;j++){
			printf("%f+j%f ,",realArr[(nLength-1)*mLength+j],imageArr[(nLength-1)*mLength+j]);
		}
		printf("%f+j%f",realArr[(nLength-1)*mLength+mLength-1],imageArr[(nLength-1)*mLength+mLength-1]);

		printf("]");
		printf("])\n");
	}
}

// power
void __vcpower1(float r,float t,float *vArr,int length,float *realArr1,float *imageArr1){

	for(int i=0;i<length;i++){
		float _r=0;

		_r=powf(r, vArr[i]);
		realArr1[i]=_r*cosf(t*vArr[i]);
		imageArr1[i]=_r*sinf(t*vArr[i]);
	}
}

/***
	b[0]*e^(j(M-1)w)+...+b[1]*e^(jw)+b[0];
	wArr1 角频率刻度
	vArr2 b/a系数
****/
void __vcpolyval(float *wArr1,int length1,float *vArr2,int length2,float *realArr3,float *imageArr3){

	for(int i=0;i<length1;i++){
		double _value1=0;
		double _value2=0;
		
		for(int j=0;j<length2;j++){ 
			_value1+=cos(wArr1[i]*(length2-1-j))*vArr2[j];
			_value2+=sin(wArr1[i]*(length2-1-j))*vArr2[j];
		}

		realArr3[i]=_value1;
		imageArr3[i]=_value2;
	}
}

void __complexMul(float real1,float image1,
				float real2,float image2,
				float *real3,float *image3){
	
	*real3=real1*real2-image1*image2;
	*image3=image1*real2+real1*image2;
}

void __complexDiv(float real1,float image1,
				float real2,float image2,
				float *real3,float *image3){
	float value=0;

	value=real2*real2+image2*image2;

	*real3=(real1*real2+image1*image2)/value;
	*image3=(image1*real2-real1*image2)/value;
}

float __complexMulM(float real1,float image1,
				float real2,float image2){
	float value=0;

	value=sqrtf(real1*real1+image1*image1)*sqrtf(real2*real2+image2*image2);

	return value;
}

float __complexDivM(float real1,float image1,
				float real2,float image2){
	float value=0;

	value=sqrtf(real1*real1+image1*image1)/sqrtf(real2*real2+image2*image2);

	return value;
}

float __complexPowM(float real,float image,int n){
	float value=0;

	value=sqrtf(real*real+image*image);
	value=powf(value, n);

	return value;
}


