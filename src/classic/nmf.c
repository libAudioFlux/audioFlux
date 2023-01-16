// 

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "nmf.h"

/***
	k k<n&&k<m,k*(n+m)<n*m
	maxIter 300
	type 1  0 Euclidean 1 KL 2 IS
	thresh 1e-3
	norm 0 max 1 sum/p-1 2 p-2
****/
void nmf(float *mDataArr,int nLength,int mLength,int k,
		float *wArr,float *hArr,
		int *maxIter,int *type,float *thresh,
		int *norm){
	float eps=1e-16;

	int _maxIter=300;
	int _type=1; // KL Div

	float _thresh=1e-3;

	int _norm=0;

	float *detW=NULL;
	float *detH=NULL;

	float *detD=NULL;
	float *dArr=NULL;

	float *wArr1=NULL; // pre
	float *hArr1=NULL;

	// cache
	float *dArr2=NULL; 
	float *dArr3=NULL; // one

	float *wArr2=NULL;
	float *hArr2=NULL;

	float *wArr3=NULL;
	float *hArr3=NULL;

	float *vArr=NULL;
	
	float value=1;

	if(maxIter){
		_maxIter=*maxIter;
	}

	if(type){
		_type=*type;
	}

	if(thresh){
		_thresh=*thresh;
	}

	if(norm){
		_norm=*norm;
	}

	detW=__vnew(nLength*k, NULL);
	detH=__vnew(k*mLength, NULL);

	detD=__vnew(nLength*mLength, NULL);
	dArr=__vnew(nLength*mLength, NULL);

	wArr1=__vnew(nLength*k, NULL);
	hArr1=__vnew(k*mLength, NULL);

	dArr2=__vnew(nLength*mLength, NULL);
	dArr3=__vnew(nLength*mLength, &value);

	wArr2=__vnew(nLength*k, NULL);
	hArr2=__vnew(k*mLength, NULL);

	wArr3=__vnew(nLength*k, NULL);
	hArr3=__vnew(k*mLength, NULL);

	vArr=__vnew(k, NULL);

	// norm W
	if(_norm==1||_norm==2){
		__mnorm(wArr, nLength, k, 0, 0, _norm, vArr);
	}
	else{
		__mmax(wArr, nLength, k, 0, vArr,NULL);
	}
	
	__mdiv_vector(wArr, vArr, 1, nLength, k, 0, wArr);

	for(int i=0;i<_maxIter;i++){
		int _t=0;

		float _h1=0;
		float _w1=0;

		memcpy(wArr1, wArr, sizeof(float )*nLength*k);
		memcpy(hArr1, hArr, sizeof(float )*k*mLength);

		// simular V
		__mdot(wArr,hArr,
			nLength,k,
			k,mLength,
			dArr);

		/***
			update H,W
			先更新H,后更新W,V误差更小些
			W归一化情况下
				IS>=KL>>Euc
			未归一化 loss不会下降!!!
				KL>Euc>>IS
		****/
		if(_type==0){ // "KL div"
			for(int j=0;j<nLength*mLength;j++){
				dArr2[j]=mDataArr[j]/(dArr[j]+eps);
			}

			// 1. updata H
			_t=2;
			__mdot2(wArr,dArr2,
				nLength,k,
				nLength,mLength,
				&_t,
				hArr2);

			__mdot2(wArr,dArr3,
				nLength,k,
				nLength,mLength,
				&_t,
				hArr3);

			for(int j=0;j<k*mLength;j++){
				hArr[j]=hArr[j]*hArr2[j]/(hArr3[j]+eps);
			}

			// 2. update W
			_t=1;
			__mdot2(dArr2,hArr,
				nLength,mLength,
				k,mLength,
				&_t,
				wArr2);

			__mdot2(dArr3,hArr,
				nLength,mLength,
				k,mLength,
				&_t,
				wArr3);

			for(int j=0;j<nLength*k;j++){
				wArr[j]=wArr[j]*wArr2[j]/(wArr3[j]+eps);
			}
		}
		else if(_type==1){ // "IS div"
			for(int j=0;j<nLength*mLength;j++){
				dArr2[j]=1.0/(dArr[j]*dArr[j]+eps)*mDataArr[j];
				dArr3[j]=1.0/(dArr[j]+eps);
			}

			// 1. updata H
			_t=2;
			__mdot2(wArr,dArr2,
				nLength,k,
				nLength,mLength,
				&_t,
				hArr2);

			__mdot2(wArr,dArr3,
				nLength,k,
				nLength,mLength,
				&_t,
				hArr3);

			for(int j=0;j<k*mLength;j++){
				hArr[j]=hArr[j]*hArr2[j]/(hArr3[j]+eps);
			}

			// 2. update W
			_t=1;
			__mdot2(dArr2,hArr,
				nLength,mLength,
				k,mLength,
				&_t,
				wArr2);

			__mdot2(dArr3,hArr,
				nLength,mLength,
				k,mLength,
				&_t,
				wArr3);

			for(int j=0;j<nLength*k;j++){
				wArr[j]=wArr[j]*wArr2[j]/(wArr3[j]+eps);
			}
		}
		else{ // "Euc dist"
			// 1. updata H
			_t=2;
			__mdot2(wArr,mDataArr,
				nLength,k,
				nLength,mLength,
				&_t,
				hArr2);

			__mdot2(wArr,dArr,
				nLength,k,
				nLength,mLength,
				&_t,
				hArr3);

			for(int j=0;j<k*mLength;j++){
				hArr[j]=hArr[j]*(hArr2[j]/(hArr3[j]+eps));
			}

			// 2. updata W
			_t=1;
			__mdot2(mDataArr,hArr,
				nLength,mLength,
				k,mLength,
				&_t,
				wArr2);

			__mdot2(dArr,hArr,
				nLength,mLength,
				k,mLength,
				&_t,
				wArr3);

			for(int j=0;j<nLength*k;j++){
				wArr[j]=wArr[j]*(wArr2[j]/(wArr3[j]+eps));
			}
		}

		// norm W
		if(_norm==1||_norm==2){
			__mnorm(wArr, nLength, k, 0, 0, _norm, vArr);
		}
		else{
			__mmax(wArr, nLength, k, 0, vArr,NULL);
		}

		__mdiv_vector(wArr, vArr, 1, nLength, k, 0, wArr);

		// check error
		__vsub(wArr, wArr1, nLength*k, detW);
		__vsub(hArr, hArr1, k*mLength, detH);
		
		_w1=__vnorm(detW, nLength*k);
		_h1=__vnorm(detH, k*mLength);
		if(_w1<_thresh&&_h1<_thresh){
			break;
		}

		// debug
		// {
		// 	float err=0;

		// 	__mdot(wArr,hArr,
		// 		nLength,k,
		// 		k,mLength,
		// 		dArr);
		// 	__vsub(mDataArr, dArr, nLength*mLength, detD);
		// 	err=__vnorm(detD, nLength*mLength);

		// 	printf("epoch %d --->	V_error: %f, W_error: %f, H_error: %f\n",
		// 			i,err,_w1,_h1);
		// }
	}

	// if(_norm){
	// 	__mmax(wArr, nLength, k, 0, vArr, NULL);

	// 	__mdiv_vector(wArr, vArr, 1, nLength, k, 0, wArr);
	// 	__mmul_vector(hArr, vArr, 0, k, mLength, 1, hArr);
	// }

	free(detW);
	free(detH);

	free(detD);
	free(dArr);

	free(wArr1);
	free(hArr1);

	free(dArr2);
	free(dArr3);

	free(wArr2);
	free(hArr2);

	free(wArr3);
	free(hArr3);

	free(vArr);
}









