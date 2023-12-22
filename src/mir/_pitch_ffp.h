

#ifndef _PITCH_FFP_H
#define _PITCH_FFP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../flux_base.h"

typedef struct OpaquePitchFFP *PitchFFPObj;

/***
	samplate 32000
	radix2Exp 12
	WindowType hamm
	slideLength (1<<radix2Exp)/4
	isContinue 0
****/
int pitchFFPObj_new(PitchFFPObj *pitchFFPObj,
					int *samplate,float *lowFre,float *highFre,
					int *radix2Exp,int *slideLength,WindowType *windowType,
					int *isContinue);

int pitchFFPObj_calTimeLength(PitchFFPObj pitchFFPObj,int dataLength);
// tempBase -18 -18/-24/...
void pitchFFPObj_setTempBase(PitchFFPObj pitchFFPObj,float tempBase);

void pitchFFPObj_pitch(PitchFFPObj pitchFFPObj,float *dataArr,int dataLength,
					float *freArr,float *dbArr);

int pitchFFPObj_getCorrData(PitchFFPObj pitchFFPObj,float **mCorrArr,float **mDbArr,float **mHeightArr,int **lenArr);
int pitchFFPObj_getCutData(PitchFFPObj pitchFFPObj,float **mCorrArr,float **mDbArr,float **mHeightArr,int **lenArr);
int pitchFFPObj_getFlagData(PitchFFPObj pitchFFPObj,int **flagArr);
int pitchFFPObj_getLightData(PitchFFPObj pitchFFPObj,float **lightArr);

int pitchFFPObj_getFormatData(PitchFFPObj pitchFFPObj,
							int **formatFlagArr,
							float **freArr1,float **freArr2,float **freArr3,
							float **dbArr1,float **dbArr2,float **dbArr3);

int pitchFFPObj_getTemporalData(PitchFFPObj pitchFFPObj,float **avgTempArr,float **maxTempArr,float **percentTempArr);

void pitchFFPObj_enableDebug(PitchFFPObj pitchFFPObj,int isDebug);
void pitchFFPObj_free(PitchFFPObj pitchFFPObj);

#ifdef __cplusplus
}
#endif

#endif