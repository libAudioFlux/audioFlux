// clang 

#ifndef DWT_FILTERCOEF_H
#define DWT_FILTERCOEF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

/***
	haar =db1
	db 2~10/20/30/40
	sym 2~10/20
	coif 1~5
	fk4,6,8,14,18,22
	bior 1.1~1.5/2.2~2.8/3.1~3.9/4.4/5.5/6.8
	coefType 0 dec 1 rec
****/
int dwt_filterCoef(WaveletDiscreteType waveletType,int t1,int t2,int coefType,
				float **loArr,float **hiArr);



#ifdef __cplusplus
}
#endif

#endif