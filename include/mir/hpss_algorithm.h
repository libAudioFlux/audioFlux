

#ifndef HPSS_ALGORITHM_H
#define HPSS_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueHPSS *HPSSObj;

/***
	windowType 'hamm'
	slideLength (1<<radix2Exp)/4
	hOrder/pOreder 21/31 (radix2Exp=11) ;must odd
****/
int hpssObj_new(HPSSObj *hpssObj,
				int radix2Exp,WindowType *windowType,int *slideLength,
				int *hOrder,int *pOrder);

int hpssObj_calDataLength(HPSSObj hpssObj,int dataLength);

void hpssObj_hpss(HPSSObj hpssObj,float *dataArr,int dataLength,float *hArr,float *pArr);

void hpssObj_free(HPSSObj hpssObj);
void hpssObj_debug(HPSSObj hpssObj);

#ifdef __cplusplus
}
#endif

#endif