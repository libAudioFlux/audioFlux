

#ifndef FLUX_COMPLEX_H
#define FLUX_COMPLEX_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>


// 复数基本运算
void __complexMul(float real1,float image1,
				float real2,float image2,
				float *real3,float *image3);
void __complexDiv(float real1,float image1,
				float real2,float image2,
				float *real3,float *image3);

float __complexMulM(float real1,float image1,
				float real2,float image2);
float __complexDivM(float real1,float image1,
				float real2,float image2);
float __complexPowM(float real,float image,int n);

// n1*m1@m1*n2
void __mcdot(float *mRealArr1,float *mImageArr1,
			float *mRealArr2,float *mImageArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			float *mRealArr3,float *mImageArr3);

// n1*m1@n2*m1
void __mcdot1(float *mRealArr1,float *mImageArr1,
			float *mRealArr2,float *mImageArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			float *mRealArr3,float *mImageArr3);

int __mcdot2(float *mRealArr1,float *mImageArr1,
			float *mRealArr2,float *mImageArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			int *type,
			float *mRealArr3,float *mImageArr3);

// 复数矩阵div
void __mcdiv(float *mRealArr1,float *mImageArr1,
			float *mRealArr2,float *mImageArr2,
			int nLength,int mLength,
			float *mRealArr3,float *mImageArr3);

void __mccut(float *mRealArr1,float *mImageArr1,
			int nLength,int mLength,
			int nIndex,int nLength1,
			int mIndex,int mLength1,
			float *mRealArr3,float *mImageArr3);

// new/convert
void __vcnew(int length,float *value,float **realArr,float **imageArr);
void __vcz2p(float *realArr,float *imageArr,int length,float *rArr,float *tArr);
void __vcp2z(float *rArr,float *tArr,int length,float *realArr,float *imageArr);

// arg/abs/conj   
void __vcarg(float *realArr1,float *imageArr1,int length,float *vArr1);
void __vcabs(float *realArr1,float *imageArr1,int length,float *vArr1);
void __vcsquare(float *realArr1,float *imageArr1,int length,float *vArr1);
void __vcconj(float *realArr1,float *imageArr1,int length,float *realArr3,float *imageArr3);

void __mcabs(float *mRealArr,float *mImageArr,int nLength,int mLength,int axis,float *mArr2);
void __mcabs1(float *mRealArr,float *mImageArr,int nLength,int mLength,int axis,int cutLength,float *mArr2);
void __mcabs2(float *mRealArr,float *mImageArr,int nLength,int mLength,int mLength2,float *mArr2);

void __mcsquare(float *mRealArr,float *mImageArr,int nLength,int mLength,int axis,float *mArr2);
void __mcsquare1(float *mRealArr,float *mImageArr,int nLength,int mLength,int axis,int cutLength,float *mArr2);
void __mcsquare2(float *mRealArr,float *mImageArr,int nLength,int mLength,int mLength2,float *mArr2);

// add/sub/mul/div
void __vcmul(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *realArr3,float *imageArr3);
void __vcdiv(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *realArr3,float *imageArr3);
void __vcadd(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *realArr3,float *imageArr3);
void __vcsub(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *realArr3,float *imageArr3);

void __vcmul_value(float *realArr1,float *imageArr1,
				float rValue,float iValue,int length,
				float *realArr3,float *imageArr3);
void __vcdiv_value(float *realArr1,float *imageArr1,
				float rValue,float iValue,int length,
				float *realArr3,float *imageArr3);
void __vcadd_value(float *realArr1,float *imageArr1,
				float rValue,float iValue,int length,
				float *realArr3,float *imageArr3);
void __vcsub_value(float *realArr1,float *imageArr1,
				float rValue,float iValue,int length,
				float *realArr3,float *imageArr3);

// mul/div 快速求模算法
void __vcmulm(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *vArr3);
void __vcdivm(float *realArr1,float *imageArr1,
			float *realArr2,float *imageArr2,int length,
			float *vArr3);
void __vcmulm_value(float *realArr1,float *imageArr1,
				float rValue,float iValue,int length,
				float *realArr3,float *imageArr3);
void __vcdivm_value(float *realArr1,float *imageArr1,
				float rValue,float iValue,int length,
				float *realArr3,float *imageArr3);

// power
void __vcpower1(float r,float t,float *vArr,int length,float *realArr1,float *imageArr1);


// polyval 多项式 阶数有小到大排列
void __vcpolyval(float *wArr1,int length1,float *vArr2,int length2,float *realArr3,float *imageArr3);

// debug相关
void __vcdebug(float *realArr1,float *imageArr1,int length,int type);
void __mcdebug(float *realArr1,float *imageArr1,int nLength,int mLength,int type);


#ifdef __cplusplus
}
#endif





#endif



