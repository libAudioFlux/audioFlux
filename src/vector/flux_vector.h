

#ifndef FLUX_VECTOR_H
#define FLUX_VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

typedef struct{
	int nLength;
	int mLength;

	int *indexArr; // _v使用nLength _m nLength+mLength

} VSlice;

// univers function op
typedef float (*UniFunc)(float element);
typedef float (*UniFunc1)(float element,float value);

// dot/mul/div/add/sub sum/min/max/mean
// n1*m1@m1*n2
void __mdot(float *mArr1,float *mArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			float *mArr3);

// n1*m1@n2*m1
void __mdot1(float *mArr1,float *mArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			float *mArr3);

int __mdot2(float *mArr1,float *mArr2,
			int nLength1,int mLength1,
			int nLength2,int mLength2,
			int *type,
			float *mArr3);

void __msub(float *mArr1,float *mArr2,int nLength,int mLength,float *mArr3);

/***
	A n*m
	1*m/n*1=>0/1/ => axis 0/1 type 0/1 vArr axis
	注:一般add/sub 同轴 mul/div 非同轴
****/
void __mmul_vector(float *mArr1,float *vArr,int type,int nLength,int mLength,int axis,float *mArr3);
void __mdiv_vector(float *mArr1,float *vArr,int type,int nLength,int mLength,int axis,float *mArr3);
void __madd_vector(float *mArr1,float *vArr,int type,int nLength,int mLength,int axis,float *mArr3);
void __msub_vector(float *mArr1,float *vArr,int type,int nLength,int mLength,int axis,float *mArr3);

void __mmul_value(float *mArr1,float value,int nLength,int mLength,float *mArr3);

// 向量归一化 type 0 p 1/2/3... 1 Inf 2 -Inf
void __mnormalize(float *mArr1,int nLength,int mLength,int axis,int type,float p,float *mArr3);
void __mnorm(float *mArr1,int nLength,int mLength,int axis,int type,float p,float *vArr2);

// order 1
void __mdiff(float *mArr1,int nLength,int mLength,int axis,int *order,float *mArr3);
// step 1
void __mdiff2(float *mArr1,int nLength,int mLength,int axis,int *step,float *mArr3);

void __mmaxfilter(float *mArr1,int nLength,int mLength,int axis,int order,float *mArr3);
void __mmedianfilter(float *mArr1,int nLength,int mLength,int axis,int order,float *mArr3);

/***
	axis 0 row 1 col -1/other all
****/
void __msum(float *mArr1,int nLength,int mLength,int axis,float *vArr3);
void __mmin(float *mArr1,int nLength,int mLength,int axis,float *vArr3,int *indexArr3);
void __mmax(float *mArr1,int nLength,int mLength,int axis,float *vArr3,int *indexArr3);
void __mmean(float *mArr1,int nLength,int mLength,int axis,float *vArr3);

void __mmedian(float *mArr1,int nLength,int mLength,int axis,float *vArr3);
void __mvar(float *mArr1,int nLength,int mLength,int axis,int type,float *vArr3);
void __mstd(float *mArr1,int nLength,int mLength,int axis,int type,float *vArr3);
void __mcov(float *mArr1,float *mArr2,int nLength,int mLength,int axis,int type,float *vArr3);
void __mcorrcoef(float *mArr1,float *mArr2,int nLength,int mLength,int axis,float *vArr3);

void __mrms(float *mArr1,int nLength,int mLength,int axis,float *vArr3);
void __menergy(float *mArr1,int nLength,int mLength,int axis,float *vArr3);
void __mzcr(float *mArr1,int nLength,int mLength,int axis,float *vArr3);

void __munwrap(float *mArr1,int nLength,int mLength,int axis);

float __vdot(float *vArr1,float *vArr2,int length);

void __vmul(float *vArr1,float *vArr2,int length,float *vArr3);
void __vdiv(float *vArr1,float *vArr2,int length,float *vArr3);
void __vadd(float *vArr1,float *vArr2,int length,float *vArr3);
void __vsub(float *vArr1,float *vArr2,int length,float *vArr3);
void __vsubi(int *vArr1,int *vArr2,int length,int *vArr3);

// order 1
void __vdiff(float *vArr1,int length,int *order,float *vArr3);
// step 1
void __vdiff2(float *vArr1,int length,int *step,float *vArr3);

void __vmaxfilter(float *vArr1,int length,int order,float *vArr3);
void __vmedianfilter(float *vArr1,int length,int order,float *vArr3);

void __vmul_value(float *vArr1,float value,int length,float *vArr3);
void __vdiv_value(float *vArr1,float value,int length,float *vArr3);
void __vadd_value(float *vArr1,float value,int length,float *vArr3);
void __vsub_value(float *vArr1,float value,int length,float *vArr3);

// 模/2范式
float __vnorm(float *vArr1,int length);
// 向量归一化 type 0 p 1/2/3... 1 Inf 2 -Inf
void __vnormalize(float *vArr1,int length,int type,float p,float *vArr2);

// 各种scale
void __vminmaxscale(float *vArr1,int length,float *vArr2);
void __vstandscale(float *vArr1,int length,int type,float *vArr2);
void __vmaxabsscale(float *vArr1,int length,float *vArr2);
void __vrobustscale(float *vArr1,int length,float *vArr2);

void __vcenterscale(float *vArr1,int length,float *vArr2);
void __vmeanscale(float *vArr1,int length,float *vArr2);
void __varctanscale(float *vArr1,int length,float *vArr2);

float __vsum(float *vArr1,int length);
int __vsumi(int *vArr1,int length);
int __vmin(float *vArr1,int length,float *value);
int __vmax(float *vArr1,int length,float *value);
int __vmini(int *vArr1,int length,int *value);
int __vmaxi(int *vArr1,int length,int *value);
int __vminabs(float *vArr1,int length,float *value);
int __vmaxabs(float *vArr1,int length,float *value);

float __vmean(float *vArr1,int length);
float __vmedian(float *vArr1,int length);
// 0 默认 N-1 1 N
float __vvar(float *vArr1,int length,int type);
float __vstd(float *vArr1,int length,int type);
float __vcov(float *vArr1,float *vArr2,int length,int type);
float __vcorrcoef(float *vArr1,float *vArr2,int length);

float __vrms(float *vArr1,int length);
float __venergy(float *vArr1,int length);
float __vzcr(float *vArr1,int length);

void __vunwrap(float *vArr1,int length,float *vArr2);

// type 0 asc 1 desc
void __vsort(float *vArr1,int length,int type,float *vArr2);
void __vsorti(int *vArr1,int length,int type,int *vArr2);

void __vcorrsort(float *vArr1,float *vArr2,float *vArr3,int length,int type);
void __vcorrsort1(float *vArr1,float *vArr2,float *vArr3,int *vArr4,int length,int type);

int __vindex(float *vArr1,int length,float value);
int __vindexi(int *vArr1,int length,int value);

int __vhas(float *vArr1,int length,float value);
int __vhasi(int *vArr1,int length,int value);

// univers function op
void __vmap(float *vArr1,int length,void *callback,float *vArr2);
void __vmap1(float *vArr1,int length,void *callback1,float value,float *vArr2);

// new/arange/linspace transpose/index相关
float *__vnew(int length,float *value);
int *__vnewi(int length,int *value);

float *__vlinspace(float start,float stop,int length,int type);
float *__vlogspace(float start,float stop,int length,int type);
int __varange(float start,float stop,float step,float **outArr);
int __varangei(int start,int stop,int step,int **outArr);

void __vfill(float *arr,int length,float value);
void __vfilli(int *arr,int length,int value);

void __vcopy(float *dstArr,float *srcArr,int length);
void __vcopyi(int *dstArr,int *srcArr,int length);

void __mtrans(float *mArr1,int nLength,int mLength,float *mArr3);

float *__vget(float *vArr1,int length,VSlice *slice);
void __vset(float *vArr1,int length,VSlice *slice,float *vArr2);

float *__mget(float *mArr1,int nLength,int mLength,VSlice *slice);
void __mset(float *mArr1,int nLength,int mLength,VSlice *slice,float *mArr2);
void __mset_value(float *mArr1,int nLength,int mLength,VSlice *slice,float value);
void __mset_vector(float *mArr1,int nLength,int mLength,VSlice *slice,int axis,float *vArr1);

// astype
void __vf2i(float *vArr1,int length,int *vArr2);
void __vi2f(int *vArr1,int length,float *vArr2);

// repeat/concat
void __vrepeat(float *vArr1,int length,int num,float *mArr3);
void __mrepeat(float *mArr1,int nLength,int mLength,int axis,int num,float *mArr3);

void __vconcat(float *vArr1,float *vArr2,int length1,int length2,float *vArr3);
void __mconcat(float *mArr1,float *mArr2,
				int nLength1,int mLength1,
				int nLength2,int mLength2,
				int axis,
				float *mArr3);

// 矩阵nLength*mLength内的裁切
void __mcut(float *mArr1,int nLength,int mLength,
			int nIndex,int nLength1,
			int mIndex,int mLength1,
			float *mArr3);

// debug vector/matrix相关 type 0 标准 1 非标准
void __vdebug(float *vArr1,int length,int type);
void __mdebug(float *mArr1,int nLength,int mLength,int type);

void __vdebugi(int *vArr1,int length,int type);
void __mdebugi(int *mArr1,int nLength,int mLength,int type);


#ifdef __cplusplus
}
#endif





#endif



