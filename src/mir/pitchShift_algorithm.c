// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "../util/flux_util.h"
#include "../dsp/resample_algorithm.h"

#include "timeStretch_algorithm.h"
#include "pitchShift_algorithm.h"

struct OpaquePitchShift{
	TimeStretchObj timeStretchObj;
	ResampleObj resampleObj;

	

};

/***
	radix2Exp 12
	WindowType hann
	slideLength (1<<radix2Exp)/4
****/
int pitchShiftObj_new(PitchShiftObj *pitchShiftObj,int *radix2Exp,int *slideLength,WindowType *windowType){
	int status=0;

	TimeStretchObj timeStretchObj=NULL;
	ResampleObj resampleObj=NULL;

	PitchShiftObj ps=NULL;

	int scale=1;
	ResampleQualityType qualType=ResampleQuality_Fast;

	ps=*pitchShiftObj=(PitchShiftObj )calloc(1,sizeof(struct OpaquePitchShift ));

	resampleObj_new(&resampleObj,&qualType,&scale,NULL);
	timeStretchObj_new(&timeStretchObj,radix2Exp,slideLength,windowType);

	ps->timeStretchObj=timeStretchObj;
	ps->resampleObj=resampleObj;

	return status;
}

// nSemitone -12~12
void pitchShiftObj_pitchShift(PitchShiftObj pitchShiftObj,int samplate,int nSemitone,float *dataArr1,int dataLength1,float *dataArr2){
	float rate=0;
	int samplate2=0;

	float *_dataArr=NULL;
	int capacity=0;
	int dataLength2=0;

	if(nSemitone>12||nSemitone<-12){
		return ;
	}

	rate=powf(2, -nSemitone*1.0/12);
	samplate2=roundf(samplate/rate);

	capacity=timeStretchObj_calDataCapacity(pitchShiftObj->timeStretchObj,rate,dataLength1);
	_dataArr=__vnew(capacity, NULL);

	// 1. timeStretch
	dataLength2=timeStretchObj_timeStretch(pitchShiftObj->timeStretchObj,rate,dataArr1,dataLength1,_dataArr);

	// 2. resample
	// resampleObj_setSamplate(pitchShiftObj->resampleObj,samplate2,samplate);
	resampleObj_setSamplateRatio(pitchShiftObj->resampleObj,rate);
	resampleObj_resample(pitchShiftObj->resampleObj,_dataArr,dataLength2,dataArr2);


	free(_dataArr);
}

void pitchShiftObj__free(PitchShiftObj pitchShiftObj){

	if(pitchShiftObj){
		timeStretchObj_free(pitchShiftObj->timeStretchObj);
		resampleObj_free(pitchShiftObj->resampleObj);

		free(pitchShiftObj);
	}
}










