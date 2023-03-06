

#ifndef REASSIGN_ALGORITHM_H
#define REASSIGN_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

typedef enum{
	Reassign_All=0,

	Reassign_Fre, 
	Reassign_Time, 

	Reassign_None, 
	
} ReassignType;

typedef struct OpaqueReassign *ReassignObj;

/***
	radix2Exp 12
	samplate 32000
	windowType 'hann'
	slideLength fftLength/4

	reType Reassign_All
	thresh 0.001 
	
	isPadding 0
	isContinue 0
****/
int reassignObj_new(ReassignObj *reassignObj,int radix2Exp,
					int *samplate,WindowType *windowType,int *slideLength,
					ReassignType *reType,float *thresh,
					int *isPadding,int *isContinue);

int reassignObj_calTimeLength(ReassignObj reassignObj,int dataLength);

// 0 complex,1 real ->mRealArr1(amp)
void reassignObj_setResultType(ReassignObj reassignObj,int type);
// order >=1
void reassignObj_setOrder(ReassignObj reassignObj,int order);

// mRealArr2&mImageArr2 may NULL
void reassignObj_reassign(ReassignObj reassignObj,float *dataArr,int dataLength,
						float *mRealArr1,float *mImageArr1,
						float *mRealArr2,float *mImageArr2);

void reassignObj_free(ReassignObj reassignObj);


#ifdef __cplusplus
}
#endif

#endif