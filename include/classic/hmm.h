

#ifndef HMM_H
#define HMM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

typedef struct OpaqueHMM *HMMObj;

int hmmObj_new(HMMObj *hmmObj,int nLength,int mLength);

void hmmObj_init(HMMObj hmmObj,float *piArr,float *mArr,float *mBArr);

// predict/decode
float hmmObj_predict(HMMObj hmmObj,int *oArr,int tLength);
float hmmObj_decode(HMMObj hmmObj,int *oArr,int tLength,int *sArr,float *mProbArr);

// maxIter 100 error 1e-3 train/genrate
void hmmObj_train(HMMObj hmmObj,int *oArr,int tLength,int *maxIter,float *error);
void hmmObj_generate(HMMObj hmmObj,int tLength,int *oArr,int *sArr);

void hmmObj_enableDebug(HMMObj hmmObj,int isDebug);

void hmmObj_free(HMMObj hmmObj);

#ifdef __cplusplus
}
#endif

#endif