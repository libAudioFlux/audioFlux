

#ifndef AUDITORY_FILTERBANK_H
#define AUDITORY_FILTERBANK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../flux_base.h"

/***
	scale liespace/mel/bark/erb/log/logspace not include edge
	scale linear include edge
	style gammatone include edge
****/
void auditory_filterBank(int num,int fftLength,int samplate,int isPseudo,
						SpectralFilterBankScaleType scaleType,
						SpectralFilterBankStyleType styleType,
						SpectralFilterBankNormalType normType,
						float lowFre,float highFre,int binPerOctave,
						float *mFilterBankArr,
						float *freBandArr,
						int *binBandArr);

// linear scale
float auditory_freToLinear(float fre,float detFre);
float auditory_linearToFre(float value,float detFre);

// linspace scale
float auditory_freToLinspace(float fre);
float auditory_linspaceToFre(float value);

// mel scale
float auditory_freToMel(float fre);
float auditory_melToFre(float mel);

// bark scale 
float auditory_freToBark(float fre);
float auditory_barkToFre(float bark);

// erb scale -> equal rect bandwidth
float auditory_freToErb(float fre);
float auditory_erbToFre(float erb);

// midi sacle
float auditory_freToMidi(float fre);
float auditory_midiToFre(float midi);

// log sacle
float auditory_freToLog(float fre,float binPerOctave);
float auditory_logToFre(float value,float binPerOctave);

// logspace scale
float auditory_freToLogspace(float fre);
float auditory_logspaceToFre(float value);

void auditory_reviseLogFre(int num,float lowFre,float highFre,int binPerOctave,int isEdge,float *lowFre3,float *highFre3);
void auditory_reviseLinearFre(int num,float lowFre,float highFre,float detFre,int isEdge,float *lowFre3,float *highFre3);

void auditory_reviseLinspaceFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3);
void auditory_reviseLogspaceFre(int num,float lowFre,float highFre,int isEdge,float *lowFre3,float *highFre3);

// length -> 4*6 matrix
float **auditory_calGammatoneCoefficient(float *freBandArr,int length,int samplate);


#ifdef __cplusplus
}
#endif

#endif