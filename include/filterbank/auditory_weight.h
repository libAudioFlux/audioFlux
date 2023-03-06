

#ifndef AUDITORY_WEIGHT_H
#define AUDITORY_WEIGHT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

void auditory_weightA(float *freArr,int length,float *dBArr);
void auditory_weightB(float *freArr,int length,float *dBArr);
void auditory_weightC(float *freArr,int length,float *dBArr);
void auditory_weightD(float *freArr,int length,float *dBArr);

#ifdef __cplusplus
}
#endif

#endif