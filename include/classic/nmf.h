

#ifndef NMF_H
#define NMF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

/***
	V=W*H
	k k<n&&k<m,k*(n+m)<n*m
	maxIter 300
	type 0 0 KL 1 IS 2 Euc 
	thresh 1e-3
	norm 0 max 1 sum/p-1 2 p-2
****/
void nmf(float *mDataArr,int nLength,int mLength,int k,
		float *wArr,float *hArr,
		int *maxIter,int *type,float *thresh,
		int *norm);


#ifdef __cplusplus
}
#endif

#endif