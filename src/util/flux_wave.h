

#ifndef FLUX_WAVE_H
#define FLUX_WAVE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

typedef struct OpaqueWaveRead *WaveReadObj;
typedef struct OpaqueWaveWrite *WaveWriteObj;

int waveReadObj_new(WaveReadObj *waveObj,char *fileName);
int waveReadObj_getInfor(WaveReadObj waveObj,int *samplate,int *bit,int *channelNum);
int waveReadObj_read(WaveReadObj waveObj,float *dataArr,int dataLength);
void waveReadObj_free(WaveReadObj waveObj);

/***
	samplate 32000
	bit 16
	channelNum 1
****/
int waveWriteObj_new(WaveWriteObj *waveObj,char *fileName,
					int *samplate,int *bit,int *channelNum);
int waveWriteObj_write(WaveWriteObj waveObj,float *dataArr,int dataLength);
void waveWriteObj_free(WaveWriteObj waveObj);

#ifdef __cplusplus
}
#endif

#endif