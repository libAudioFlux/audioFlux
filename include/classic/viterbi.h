

#ifndef VITERBI_H
#define VITERBI_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

/***
	λ=(π,A,B) sLength*1,sLength*sLength,sLength*nLength
	oArr B nLength index
	mProbArr tLength*sLength
	sArr A sLength index
****/
float viterbi(float *piArr,float *mAArr,float *mBArr,
			int sLength,int nLength,
			int *oArr,int tLength,int *isLog,
			int *sArr,float *mProbArr,int *mIndexArr);


#ifdef __cplusplus
}
#endif

#endif