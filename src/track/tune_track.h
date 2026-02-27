

#ifndef TUNE_TRACK_H
#define TUNE_TRACK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

typedef struct OpaqueTuneTrack *TuneTrackObj;

/***
	samplate 32000
	lowFre 27
	highFre 4000
	radix2Exp 12
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int tuneTrackObj_new(TuneTrackObj *tuneTrackObj,
					int *samplate,float *lowFre,float *highFre,
					int *radix2Exp,int *slideLength,int *isContinue);

int tuneTrackObj_calTimeLength(TuneTrackObj tuneTrackObj,int dataLength);
void tuneTrackObj_clear(TuneTrackObj tuneTrackObj);

// tempBase -18 -18/-24/...
void tuneTrackObj_setTempBase(TuneTrackObj tuneTrackObj,float tempBase);
// 5/8 ->220, base>=1
void tuneTrackObj_setUpdateBase(TuneTrackObj tuneTrackObj,float minBase,float maxBase);

void tuneTrackObj_tune(TuneTrackObj tuneTrackObj,float *dataArr,int dataLength,
					float *freArr);

int tuneTrackObj_getDataArr(TuneTrackObj tuneTrackObj,float **valueArr,float **dbArr,int **countArr,int **flagArr,float **fluxArr);

int tuneTrackObj_getTemporalDataArr(TuneTrackObj tuneTrackObj,float **avgArr,float **maxArr,float **percentArr);

void tuneTrackObj_enableDebug(TuneTrackObj tuneTrackObj,int isDebug);
void tuneTrackObj_free(TuneTrackObj tuneTrackObj);

#ifdef __cplusplus
}
#endif

#endif