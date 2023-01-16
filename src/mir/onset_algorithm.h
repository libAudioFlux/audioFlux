

#ifndef ONSET_ALGORITHM_H
#define ONSET_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

typedef enum{
	Novelty_Flux=0,

	Novelty_HFC,
	Novelty_SD,
	Novelty_SF,
	Novelty_MKL,

	Novelty_PD,
	Novelty_WPD, 
	Novelty_NWPD, 

	Novelty_CD, 
	Novelty_RCD, 

	Novelty_Broadband,

} NoveltyType;

typedef struct {
	int step; // >=1
	float p; // 1 !=0
	int isPostive; // 1
	int isExp; // 0
	int type; // 0 sum 1 mean

	float threshold; // >=0

	int isNorm; // 0|1
	float gamma; // 1 0.5/1/10/20...

} NoveltyParam;

typedef struct OpaqueOnset *OnsetObj;

/***
	samplate 32000
	filterOrder 1
	type Novelty_Flux
****/
int onsetObj_new(OnsetObj *onsetObj,int nLength,int mLength,int slideLength,
				int *samplate,int *filterOrder,
				NoveltyType *type);

int onsetObj_onset(OnsetObj onsetObj,float *mDataArr1,float *mDataArr2,
				NoveltyParam *param,int *indexArr,int indexLength,
				float *evnArr,int *pointArr);

void onsetObj_free(OnsetObj onsetObj);
void onsetObj_debug(OnsetObj onsetObj);

#ifdef __cplusplus
}
#endif

#endif