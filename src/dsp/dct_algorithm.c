// clang -g -c 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"

#include "dct_algorithm.h"

typedef enum{
	DCTExec_None=0,

	DCTExec_DCT,
	DCTExec_IDCT,

} DCTExecType;

struct OpaqueDCT{
	DCTExecType execType; 

	int length; // 数据长度

	float s0,s1; // scale系数
	float *mArr; // cos系数 k*n

	float *resultArr; 
	float *normalArr;
};

// getResultData->copyResultData 废弃掉
// int dctObj_getResultData(DCTObj dctObj,float **dataArr,float **normalArr);
// int dctObj_copyResultData(DCTObj dctObj,float **dataArr,float **normalArr);
// void dctObj_debug(DCTObj dctObj);

static float *__calCosMatrix(int length);
static void _vdot(float **arrArr,float *dataArr,float *outArr,int length);

int dctObj_new(DCTObj *dctObj,int length,DCTType *type){
	int status=0;
	DCTObj dct=NULL;

	float s0=0,s1=0;

	float *mArr=NULL;
	float *resultArr=NULL;
	float *normalArr=NULL;

	if(length<=0){
		return -100;
	}

	dct=*dctObj=(DCTObj )calloc(1, sizeof(struct OpaqueDCT ));

	s0=sqrtf(1.0/length);
	s1=sqrtf(2.0/length);

	resultArr=(float *)calloc(length, sizeof(float ));
	normalArr=(float *)calloc(length, sizeof(float ));
	mArr=__calCosMatrix(length);

	// printf("DCTMatrix is:\n");
	// __mdebug(mArr, length, length, 1);

	dct->length=length;

	dct->s0=s0;
	dct->s1=s1;

	dct->resultArr=resultArr;
	dct->normalArr=normalArr;
	dct->mArr=mArr;

	return status;
}

/***
	DCT-II 朴素实现算法
	X(k)=sqrt(2/N)s(k)∑x(n)cos(PI/N*(n+0.5)*k)
	k=0,s(k)=sqrt(0.5);other s(k)=1
****/
void dctObj_dct(DCTObj dctObj,float *dataArr1,int isNorm,float *dataArr2){
	int length=0;

	float s0=0,s1=0;

	float *mArr=NULL;
	float *resultArr=NULL;
	float *normalArr=NULL;

	length=dctObj->length;
	s0=dctObj->s0;
	s1=dctObj->s1;

	mArr=dctObj->mArr;
	resultArr=dctObj->resultArr;
	normalArr=dctObj->normalArr;

	__mdot(mArr,dataArr1,
		length,length,
		length,1,
		dataArr2);

	if(isNorm){
		dataArr2[0]=dataArr2[0]*s0;
		for(int i=1;i<length;i++){
			dataArr2[i]=dataArr2[i]*s1;
		}
	}
}

void dctObj_idct(DCTObj dctObj,float *dataArr1,int isNorm,float *dataArr2){

}

int dctObj_getResultData(DCTObj dctObj,float **dataArr,float **normalArr){
	int length=0;

	length=dctObj->length;
	*dataArr=dctObj->resultArr;
	*normalArr=dctObj->normalArr;

	return length;
}

int dctObj_copyResultData(DCTObj dctObj,float **dataArr,float **normalArr){
	int length=0;

	float *dArr=NULL;
	float *nArr=NULL;

	length=dctObj->length;
	dArr=(float *)calloc(length, sizeof(float ));
	nArr=(float *)calloc(length, sizeof(float ));
	for(int i=0;i<length;i++){
		dArr[i]=dctObj->resultArr[i];
		nArr[i]=dctObj->normalArr[i];
	}

	*dataArr=dArr;
	*normalArr=nArr;

	return length;
}

void dctObj_free(DCTObj dctObj){
	int length=0;
	float *mArr=NULL; // cos系数

	float *resultArr=NULL; 
	float *normalArr=NULL;

	if(!dctObj){
		return;
	}

	length=dctObj->length;
	mArr=dctObj->mArr;

	resultArr=dctObj->resultArr;
	normalArr=dctObj->normalArr;

	free(mArr);

	free(resultArr);
	free(normalArr);

	free(dctObj);
}

static float *__calCosMatrix(int length){
	float *mArr=NULL;

	mArr=__vnew(length*length, NULL);
	for(int i=0;i<length;i++){
		for(int j=0;j<length;j++){
			mArr[i*length+j]=cosf(M_PI*(j+0.5)*i/length);
		}
	}

	return mArr;
}

static void _vdot(float **arrArr,float *dataArr,float *outArr,int length){
	float value=0;

	for(int i=0;i<length;i++){
		value=0;
		for(int j=0;j<length;j++){
			value+=arrArr[i][j]*dataArr[j];
		}
		outArr[i]=value;
	}
}

void dctObj_debug(DCTObj dctObj){
	int length=0;

	float *resultArr=NULL; 
	float *normalArr=NULL;

	length=dctObj->length;

	resultArr=dctObj->resultArr;
	normalArr=dctObj->normalArr;

	printf("dct result is:\n");
	for(int i=0;i<length;i++){
		printf("%f,",resultArr[i]);
	}
	printf("\n");

	printf("dct normal is:\n");
	for(int i=0;i<length;i++){
		printf("%f,",normalArr[i]);
	}
	printf("\n");
}





