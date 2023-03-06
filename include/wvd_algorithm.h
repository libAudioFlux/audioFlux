

#ifndef WVD_ALGORITHM_H
#define WVD_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueWVD *WVDObj;

int wvdObj_new(WVDObj *wvdObj,int num,int radix2Exp);

void wvdObj_wvd(WVDObj wvdObj,float *dataArr,float *mRealArr1,float *mImageArr1);

void wvdObj_free(WVDObj wvdObj);

#ifdef __cplusplus
}
#endif

#endif