

#ifndef SPECTRAL_ALGORITHM_H
#define SPECTRAL_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

typedef struct OpaqueSpectral *SpectralObj;

int spectralObj_new(SpectralObj *spectralObj,int num,float *freBandArr);

// 0~num-1
void spectralObj_setEdge(SpectralObj spectralObj,int start,int end);
void spectralObj_setEdgeArr(SpectralObj spectralObj,int *indexArr,int indexLength);

void spectralObj_setTimeLength(SpectralObj spectralObj,int timeLength);

void spectralObj_flatness(SpectralObj spectralObj,float *mDataArr,float *dataArr);
// isExp 0;step 1,p 1,type 0 sum 1 mean;
void spectralObj_flux(SpectralObj spectralObj,float *mDataArr,
					int step,float p,int isPostive,int *isExp,int *type,
					float *dataArr);
void spectralObj_rolloff(SpectralObj spectralObj,float *mDataArr,float threshold,float *dataArr);

void spectralObj_centroid(SpectralObj spectralObj,float *mDataArr,float *dataArr);
void spectralObj_spread(SpectralObj spectralObj,float *mDataArr,float *dataArr);
void spectralObj_skewness(SpectralObj spectralObj,float *mDataArr,float *dataArr);
void spectralObj_kurtosis(SpectralObj spectralObj,float *mDataArr,float *dataArr);

void spectralObj_entropy(SpectralObj spectralObj,float *mDataArr,int isNorm,float *dataArr);
void spectralObj_crest(SpectralObj spectralObj,float *mDataArr,float *dataArr);
void spectralObj_slope(SpectralObj spectralObj,float *mDataArr,float *dataArr);
void spectralObj_decrease(SpectralObj spectralObj,float *mDataArr,float *dataArr);
// p 2 !=0
void spectralObj_bandWidth(SpectralObj spectralObj,float *mDataArr,float p,float *dataArr);
void spectralObj_rms(SpectralObj spectralObj,float *mDataArr,float *dataArr);
// gamma 10
void spectralObj_energy(SpectralObj spectralObj,float *mDataArr,int isLog,float gamma,float *dataArr);

void spectralObj_hfc(SpectralObj spectralObj,float *mDataArr,float *dataArr);
void spectralObj_sd(SpectralObj spectralObj,float *mDataArr,int step,int isPostive,float *dataArr);
void spectralObj_sf(SpectralObj spectralObj,float *mDataArr,int step,int isPostive,float *dataArr);
void spectralObj_mkl(SpectralObj spectralObj,float *mDataArr,int type,float *dataArr);

void spectralObj_pd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr);
void spectralObj_wpd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr);
void spectralObj_nwpd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr);

void spectralObj_cd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr);
void spectralObj_rcd(SpectralObj spectralObj,float *mSpecArr,float *mPhaseArr,float *dataArr);

// mag|power threshold >=0
void spectralObj_broadband(SpectralObj spectralObj,float *mDataArr,float threshold,float *dataArr);

/***
	step >=1
	threshold 0
	methodType 'sub'
	dataType 'value'
****/
void spectralObj_novelty(SpectralObj spectralObj,float *mDataArr,
						int step,float threshold,
						SpectralNoveltyMethodType *methodType,SpectralNoveltyDataType *dataType,
						float *dataArr);

void spectralObj_eef(SpectralObj spectralObj,float *mDataArr,int isNorm,float *dataArr);
// gamma 1/10/20... song 0.5
void spectralObj_eer(SpectralObj spectralObj,float *mDataArr,int isNorm,float gamma,float *dataArr);

// statistics
void spectralObj_max(SpectralObj spectralObj,float *mDataArr,float *valueArr,float *freArr);
void spectralObj_mean(SpectralObj spectralObj,float *mDataArr,float *valueArr,float *freArr);
void spectralObj_var(SpectralObj spectralObj,float *mDataArr,float *valueArr,float *freArr);

void spectralObj_free(SpectralObj spectralObj);

#ifdef __cplusplus
}
#endif

#endif