

#ifndef CHROMA_ALGORITHM_H
#define CHROMA_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueChroma *ChromaObj;

// chNum >=binPerOctave; times
int chromaObj_new(ChromaObj *chromaObj,int num,float *freBandArr,
				int *binPerOctave,
				int *chNum);

void chromaObj_setTimeLength(ChromaObj chromaObj,int timeLength);

// mDataArr1 ->timeLength*num; mDataArr2 ->timeLength*chNum
void chromaObj_chroma(ChromaObj chromaObj,float *mDataArr1,float *mDataArr2);

void chromaObj_free(ChromaObj chromaObj);

#ifdef __cplusplus
}
#endif

#endif