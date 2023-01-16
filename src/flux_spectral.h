

#ifndef FLUX_SPECTRAL_H
#define FLUX_SPECTRAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flux_base.h"

void spectral_flatness(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,
					float *vArr);

// step >=1;type 0(default) sum 1 mean
void spectral_flux(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					int step,float p,int isPostive,int isExp,int type,
					float *vArr);

void spectral_rolloff(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,float threshold,
					float *vArr);

void spectral_centroid(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,
					float *vArr);

void spectral_spread(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,float *cArr,
					float *vArr);

void spectral_skewness(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,float *cArr1,float *cArr2,
					float *vArr);

void spectral_kurtosis(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,float *cArr1,float *cArr2,
					float *vArr);

// 1 matlab 0 song
void spectral_entropy(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *sumArr,int isNorm,
					float *vArr);

void spectral_crest(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *sumArr,
					float *vArr);

void spectral_slope(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *meanFreArr,float *meanValueArr,
					float *vArr);

void spectral_decrease(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *sumArr,
					float *vArr);

void spectral_bandWidth(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float *freArr,float *cArr,float p,
					float *vArr);

void spectral_rms(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr);


// ['hfc', 'sd', 'sf', 'mkl', 'pd', 'wpd', 'nwpd', 'cd', 'rcd']
void spectral_hfc(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr);

void spectral_sd(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				int step,int isPostive,
				float *vArr);

void spectral_sf(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				int step,int isPostive,
				float *vArr);

// type 0 sum 1 mean
void spectral_mkl(float *mDataArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				int type,
				float *vArr);

void spectral_pd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr);

void spectral_wpd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr);

void spectral_nwpd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr);

void spectral_cd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr);

void spectral_rcd(float *mSpecArr,float *mPhaseArr,int nLength,int mLength,
				int *indexArr,int indexLength,
				float *vArr);

// threshold 0/3/6  
void spectral_broadband(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					float threshold,
					float *vArr);

/***
	step >=1
	threshold 0
	methodType 'sub'
	dataType 'value'
****/
void spectral_novelty(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					int step,float threshold,
					SpectralNoveltyMethodType *methodType,SpectralNoveltyDataType *dataType,
					float *vArr);

void spectral_energy(float *mDataArr,int nLength,int mLength,
					int *indexArr,int indexLength,
					int isPower,int isLog,float gamma,
					float *vArr);


#ifdef __cplusplus
}
#endif

#endif