

#ifndef NSGT_ALGORITHM_H
#define NSGT_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef enum{
	NSGTFilterBank_Efficient=0, 
	NSGTFilterBank_Standard, 

} NSGTFilterBankType;

typedef struct OpaqueNSGT *NSGTObj;

/***
	samplate 32000
	lowFre linear 0 mel/bark/erb 27.5 log 32.703
	highFre 16000
	binPerOctave 12 12/24/36
	minLen 3 >=1
	
	nsgtFilterBankType "efficient"
	filterStyleType "hann" !Gammatone
	filterScaleType "log" Linear~Log
	filterNormalType "bandWidth" None||BandWidth
****/
int nsgtObj_new(NSGTObj *nsgtObj,int num,int radix2Exp,
				int *samplate,float *lowFre,float *highFre,int *binPerOctave,
				int *minLen,
				NSGTFilterBankType *nsgtFilterBankType,
				SpectralFilterBankScaleType *filterScaleType,
				SpectralFilterBankStyleType *filterStyleType,
				SpectralFilterBankNormalType *filterNormalType);

int nsgtObj_getMaxTimeLength(NSGTObj nsgtObj);
int nsgtObj_getTotalTimeLength(NSGTObj nsgtObj);
int *nsgtObj_getTimeLengthArr(NSGTObj nsgtObj);

float *nsgtObj_getFreBandArr(NSGTObj nsgtObj);
int *nsgtObj_getBinBandArr(NSGTObj nsgtObj);

// default 3  minLength>=1
void nsgtObj_setMinLength(NSGTObj nsgtObj,int minLength);

void nsgtObj_nsgt(NSGTObj nsgtObj,float *dataArr,float *mRealArr3,float *mImageArr3);

// test cell data
void nsgtObj_getCellData(NSGTObj nsgtObj,float **realArr3,float **imageArr3);

void nsgtObj_free(NSGTObj nsgtObj);

#ifdef __cplusplus
}
#endif

#endif