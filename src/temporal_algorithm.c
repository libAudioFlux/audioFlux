// clang 

#include <string.h>
#include <math.h>

#include "vector/flux_vector.h"
#include "vector/flux_vectorOp.h"
#include "vector/flux_complex.h"

#include "util/flux_util.h"

#include "dsp/flux_window.h"
#include "dsp/fft_algorithm.h"

#include "temporal_algorithm.h"

struct OpaqueTemporal{
	int frameLength;
	int slideLength;
	
	float *windowDataArr; // frameLength

	int timeLength;
	float *mDataArr; // timeLength*frameLength

	float *energyArr; // energy x(t)^2
	float *rmsArr; // rms sqrt(E/N)
	float *zcrArr; // zero-crossRate 

};

/***
	frameLength 2048
	slideLength 512
	windowType Hann
****/
int temporalObj_new(TemporalObj *temporalObj,
			 	int *frameLength,int *slideLength,WindowType *windowType){
	int status=0;

	int _frameLength=2048;
	int _slideLength=0;

	WindowType _windowType=Window_Hann;

	TemporalObj tmp=NULL;
	float *windowDataArr=NULL;

	tmp=*temporalObj=(TemporalObj )calloc(1, sizeof(struct OpaqueTemporal ));

	if(frameLength){
		if(*frameLength>0){
			_frameLength=*frameLength;
		}
	}

	_slideLength=_frameLength/4;
	if(slideLength){
		if(*slideLength>0){
			_slideLength=*slideLength;
		}
	}

	if(windowType){
		_windowType=*windowType;
	}

	windowDataArr=window_calFFTWindow(_windowType,_frameLength);

	tmp->frameLength=_frameLength;
	tmp->slideLength=_slideLength;

	tmp->windowDataArr=windowDataArr;

	return status;
}

int temporalObj_calTimeLength(TemporalObj temporalObj,int dataLength){
	int frameLength=0; 
	int slideLength=0;
	int timeLength=0;

	frameLength=temporalObj->frameLength;
	slideLength=temporalObj->slideLength;
	if(dataLength<frameLength){
		return 0;
	}

	timeLength=(dataLength-frameLength)/slideLength+1;
	return timeLength;
}

void temporalObj_temporal(TemporalObj temporalObj,float *dataArr,int dataLength){
	int frameLength=0; 
	int slideLength=0;
	int timeLength=0;

	float *windowDataArr=NULL; // frameLength

	float *energyArr=NULL; // energy x(t)^2
	float *rmsArr=NULL; // rms sqrt(E/N)
	float *zcrArr=NULL;

	float *mDataArr=NULL;

	frameLength=temporalObj->frameLength;
	slideLength=temporalObj->slideLength;
	if(dataLength<frameLength){
		return ;
	}

	timeLength=(dataLength-frameLength)/slideLength+1;
	if(temporalObj->timeLength<timeLength||
		temporalObj->timeLength>2*timeLength){

		free(temporalObj->energyArr);
		free(temporalObj->rmsArr);
		free(temporalObj->zcrArr);

		free(temporalObj->mDataArr);

		temporalObj->energyArr=__vnew(timeLength, NULL);
		temporalObj->rmsArr=__vnew(timeLength, NULL);
		temporalObj->zcrArr=__vnew(timeLength, NULL);

		temporalObj->mDataArr=__vnew(timeLength*frameLength, NULL);
	}

	windowDataArr=temporalObj->windowDataArr;

	energyArr=temporalObj->energyArr;
	rmsArr=temporalObj->rmsArr;
	zcrArr=temporalObj->zcrArr;

	mDataArr=temporalObj->mDataArr;

	// energy/rms/zcr
	for(int i=0;i<timeLength;i++){
		__vmul(dataArr+i*slideLength, windowDataArr, frameLength, mDataArr+i*frameLength);

		energyArr[i]=__venergy(mDataArr+i*frameLength, frameLength);
		rmsArr[i]=sqrtf(energyArr[i]/frameLength);
		zcrArr[i]=__vzcr(mDataArr+i*frameLength, frameLength);
	}

	temporalObj->timeLength=timeLength;
}

void temporalObj_getData(TemporalObj temporalObj,
						float **eArr,float **rArr,float **zArr,
						float **mDataArr){
	if(eArr){
		*eArr=temporalObj->energyArr;
	}

	if(rArr){
		*rArr=temporalObj->rmsArr;
	}

	if(zArr){
		*zArr=temporalObj->zcrArr;
	}

	if(mDataArr){
		*mDataArr=temporalObj->mDataArr;
	}
}

void temporalObj_ezr(TemporalObj temporalObj,float gamma,float *vArr3){
	float *energyArr=NULL; 
	float *zcrArr=NULL;

	int frameLength=0;
	int timeLength=0;

	float value1=0;
	float value2=0;

	energyArr=temporalObj->energyArr;
	zcrArr=temporalObj->zcrArr;

	frameLength=temporalObj->frameLength;
	timeLength=temporalObj->timeLength;
	for(int i=0;i<timeLength;i++){
		value1=log10f(1+energyArr[i]*gamma);
		value2=zcrArr[i]*frameLength+1;

		vArr3[i]=value1/value2;
	}
}

void temporal_free(TemporalObj temporalObj){

	if(temporalObj){
		free(temporalObj->windowDataArr);
		free(temporalObj->mDataArr);

		free(temporalObj->energyArr);
		free(temporalObj->rmsArr);
		free(temporalObj->zcrArr);

		free(temporalObj);
	}
}









