

#ifndef _PITCH_UTIL_H
#define _PITCH_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

float __queue_fre3(float value1,float value2,float value3,
					int *s1,int *s2,
					int *k1,int *k2,int *k3);

float __queue_fre2(float value1,float value2,
					int *k1,int *k2);

#ifdef __cplusplus
}
#endif

#endif