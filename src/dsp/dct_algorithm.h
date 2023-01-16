

#ifndef DCT_ALGORITHM_H
#define DCT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef enum{
	DCT_2=0, // 默认

	DCT_1,
	DCT_3,
	DCT_4,
	DCT_5,
	DCT_6,
	DCT_7,
	DCT_8,

} DCTType;

typedef struct OpaqueDCT *DCTObj;

int dctObj_new(DCTObj *dctObj,int length,DCTType *type);

void dctObj_dct(DCTObj dctObj,float *dataArr1,int isNorm,float *dataArr2);
void dctObj_idct(DCTObj dctObj,float *dataArr1,int isNorm,float *dataArr2);

void dctObj_free(DCTObj dctObj);

#ifdef __cplusplus
}
#endif

#endif