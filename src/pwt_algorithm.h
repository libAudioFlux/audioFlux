

#ifndef PWT_ALGORITHM_H
#define PWT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef struct OpaquePWT *PWTObj;

int pwtObj_new(PWTObj *pwtObj,int num,int radix2Exp,
			 int *samplate,float *lowFre,float *highFre,int *binPerOctave,
			 SpectralFilterBankScaleType *scaleType,
			 SpectralFilterBankStyleType *styleType,
			 SpectralFilterBankNormalType *normalType,
			 int *isPadding);

float *pwtObj_getFreBandArr(PWTObj pwtObj);
int *pwtObj_getBinBandArr(PWTObj pwtObj);

void pwtObj_pwt(PWTObj pwtObj,float *dataArr,float *mRealArr3,float *mImageArr3);

void pwtObj_enableDet(PWTObj pwtObj,int flag);
void pwtObj_pwtDet(PWTObj pwtObj,float *dataArr,float *mRealArr3,float *mImageArr3);

void pwtObj_free(PWTObj pwtObj);

#ifdef __cplusplus
}
#endif

#endif