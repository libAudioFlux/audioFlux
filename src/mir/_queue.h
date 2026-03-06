

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

int __queue_num(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,int **indexArr2,int **kArr,int **numArr);
float __queue_multi(float *freArr,float *dbArr,float *heightArr,int length,int num,int subType,int unionType,int direction);

int __queue_bear(float *freArr,float *dbArr,float *heightArr,int length,float min,float base,int *index);
int __queue_count(float *freArr,float *dbArr,float *heightArr,int length,int start,float min,float base,int step);

float __queue_standard(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,
					float *freArr2,float *dbArr2,float *heightArr2,int length2,
					float *freArr3,float *dbArr3,float *heightArr3,int refLength,
					float light,int *valid,
					int *formatFlag,
					float *fre1,float *fre2,float *fre3,
					float *db1,float *db2,float *db3);

float __queue_cut(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,
				float *freArr2,float *dbArr2,float *heightArr2,int length2,
				float *freArr3,float *dbArr3,float *heightArr3,int refLength,
				float light,int *valid,
				int *formatFlag,
				float *fre1,float *fre2,float *fre3,
				float *db1,float *db2,float *db3);

float __queue_fast(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,
				float *freArr2,float *dbArr2,float *heightArr2,int refLength,
				float light,int *valid,
				int *formatFlag,
				float *fre1,float *fre2,float *fre3,
				float *db1,float *db2,float *db3);

float __queue_direct(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,float light,int *valid);
float __queue_slide(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,float light,int *valid,int *status);
float __queue_weak(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,float light,int *valid,int *status);

#ifdef __cplusplus
}
#endif

#endif