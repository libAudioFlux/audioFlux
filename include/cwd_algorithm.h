

#ifndef CWD_ALGORITHM_H
#define CWD_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaqueCWD *CWDObj;

int cwdObj_new(CWDObj *wvdObj,int num,int radix2Exp);

void cwdObj_cwd(CWDObj wvdObj,float *dataArr,float *mRealArr1,float *mImageArr1);

void cwdObj_free(CWDObj wvdObj);

#ifdef __cplusplus
}
#endif

#endif