

#ifndef _TRIST3_H
#define _TRIST3_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

int trist3(float *corrArr1,float *dbArr1,float *heightArr1,int length1,
			float *corrArr2,float *dbArr2,float *heightArr2,int length2,
			float *corrArr3,float *dbArr3,float *heightArr3,int length3,
			float light,float *outFre,
			int *formatFlag,
			float *fre1,float *fre2,float *fre3,
			float *db1,float *db2,float *db3);


#ifdef __cplusplus
}
#endif

#endif