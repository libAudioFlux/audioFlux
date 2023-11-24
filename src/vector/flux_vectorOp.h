

#ifndef FLUX_VECTOROP_H
#define FLUX_VECTOROP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

// element math相关
// abs/ng/floor/ceil/round vector==matrix
void __vabs(float *vArr1,int length,float *vArr2);
void __vng(float *vArr1,int length,float *vArr2);
void __vfloor(float *vArr1,int length,float *vArr2);
void __vceil(float *vArr1,int length,float *vArr2);
void __vround(float *vArr1,int length,float *vArr2);

// cos/sin/tan/acos/asin/atan
void __vcos(float *vArr1,int length,float *vArr2);
void __vsin(float *vArr1,int length,float *vArr2);
void __vtan(float *vArr1,int length,float *vArr2);
void __vacos(float *vArr1,int length,float *vArr2);
void __vasin(float *vArr1,int length,float *vArr2);
void __vatan(float *vArr1,int length,float *vArr2);

// exp/exp2/pow/sqrt
void __vexp(float *vArr1,int length,float *vArr2);
void __vexp2(float *vArr1,int length,float *vArr2);
void __vpow(float *vArr1,float exp,int length,float *vArr2);
void __vsqrt(float *vArr1,int length,float *vArr2);

// log/log10/log2
void __vlog(float *vArr1,int length,float *vArr2);
void __vlog10(float *vArr1,int length,float *vArr2);
void __vlog2(float *vArr1,int length,float *vArr2);

// log compress
void __vlog_compress(float *vArr1,float *gamma,float *base,int length,float *vArr2);
void __vlog10_compress(float *vArr1,float *gamma,float *base,int length,float *vArr2);
void __vlog2_compress(float *vArr1,float *gamma,float *base,int length,float *vArr2);

// sinc
void __vsinc(float *vArr1,int length,float *vArr2);
void __vsinc_low(float *vArr1,int length,float cut,float *vArr2);
void __vsinc_high(float *vArr1,int length,float cut,float *vArr2);
void __vsinc_bandpass(float *vArr1,int length,float cut1,float cut2,float *vArr2);
void __vsinc_stop(float *vArr1,int length,float cut1,float cut2,float *vArr2);

// gradient edge 1/2
void __vgradient(float *vArr1,int length,int edgeOrder,float *vArr2);

// interpolation
void __vinterp_linear(float *xArr1,float *yArr1,int length1,float *xArr2,int length2,float *yArr2);

// pad type 0 'constant' 1 'refect' 2 'wrap' 3 'symmetric' 4 'edge'
void __vpad(float *vArr1,int length,
			int headLen,int tailLen,int type,
			float *value1,float *value2,
			float *vArr2);

void __mpad(float *mArr1,int nlength,int mLength,
			int headLen,int tailLen,int type,
			float *value1,float *value2,
			int axis,
			float *mArr2);

// constant
void __vpad_center1(float *vArr1,int vLength,
				int leftLength,int rightLength,
				float value1,float value2);
void __vpad_left1(float *vArr1,int vLength,int leftLength,int value);
void __vpad_right1(float *vArr1,int vLength,int rightLength,int value);

// reflect
void __vpad_center2(float *vArr1,int vLength,int leftLength,int rightLength);
void __vpad_left2(float *vArr1,int vLength,int leftLength);
void __vpad_right2(float *vArr1,int vLength,int rightLength);

// wrap
void __vpad_center3(float *vArr1,int vLength,int leftLength,int rightLength);
void __vpad_left3(float *vArr1,int vLength,int leftLength);
void __vpad_right3(float *vArr1,int vLength,int rightLength);

// symmetric
void __vpad_center4(float *vArr1,int vLength,int leftLength,int rightLength);
void __vpad_left4(float *vArr1,int vLength,int leftLength);
void __vpad_right4(float *vArr1,int vLength,int rightLength);

// edge
void __vpad_center5(float *vArr1,int vLength,int leftLength,int rightLength);
void __vpad_left5(float *vArr1,int vLength,int leftLength);
void __vpad_right5(float *vArr1,int vLength,int rightLength);



#ifdef __cplusplus
}
#endif

#endif