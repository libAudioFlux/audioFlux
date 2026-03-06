// 

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "../util/flux_util.h"

#include "_queue.h"
#include "_trist3.h"

static int __trist_standard(float *corrArr1,float *dbArr1,float *heightArr1,int length1,
							float *corrArr2,float *dbArr2,float *heightArr2,int length2,
							float *corrArr3,float *dbArr3,float *heightArr3,int refLength,
							float light,float *outFre,int *valid,
							int *formatFlag,
							float *fre1,float *fre2,float *fre3,
							float *db1,float *db2,float *db3);

static int __trist_cut(float *corrArr1,float *dbArr1,float *heightArr1,int length1,
					float *corrArr2,float *dbArr2,float *heightArr2,int length2,
					float *corrArr3,float *dbArr3,float *heightArr3,int refLength,
					float light,float *outFre,int *valid,
					int *formatFlag,
					float *fre1,float *fre2,float *fre3,
					float *db1,float *db2,float *db3);

static int __trist_fast(float *corrArr1,float *dbArr1,float *heightArr1,int length1,
						float *corrArr2,float *dbArr2,float *heightArr2,int refLength,
						float light,float *outFre,int *valid,
						int *formatFlag,
						float *fre1,float *fre2,float *fre3,
						float *db1,float *db2,float *db3);

static int __trist(float *corrArr1,float *dbArr1,float *heightArr1,int length1,float light,float *outFre,int *valid,int *status);

static int __arr_maxIndex(float *arr,int length);

// queue_slide
int trist3(float *corrArr1,float *dbArr1,float *heightArr1,int length1,
			float *corrArr3,float *dbArr3,float *heightArr3,int length3,
			float *corrArr5,float *dbArr5,float *heightArr5,int length5,
			float light,float *outFre,
			int *formatFlag,
			float *fre1,float *fre2,float *fre3,
			float *db1,float *db2,float *db3){

	int flag=0;
	float fre=0;

	int valid=0;
	int status=0;

	if(!fre&&length5){
		flag=__trist_standard(corrArr5,dbArr5,heightArr5,length5,
							corrArr3,dbArr3,heightArr3,length3,
							corrArr1,dbArr1,heightArr1,length1,
							light,&fre,&valid,
							formatFlag,
							fre1,fre2,fre3,
							db1,db2,db3);
	}
	
	if(!fre&&length5){
		flag=__trist_cut(corrArr5,dbArr5,heightArr5,length5,
						corrArr3,dbArr3,heightArr3,length3,
						corrArr1,dbArr1,heightArr1,length1,
						light,&fre,&valid,
						formatFlag,
						fre1,fre2,fre3,
						db1,db2,db3);
	}

	if(!fre&&length3){
		flag=__trist_fast(corrArr3,dbArr3,heightArr3,length3,
						corrArr1,dbArr1,heightArr1,length1,
						light,&fre,&valid,
						formatFlag,
						fre1,fre2,fre3,
						db1,db2,db3);
	}
	
	if(!fre&&length1){
		flag=__trist(corrArr1,dbArr1,heightArr1,length1,light,&fre,&valid,&status);
	}
	
	*outFre=fre;
	return flag;
}

static int __trist_standard(float *corrArr1,float *dbArr1,float *heightArr1,int length1,
							float *corrArr3,float *dbArr3,float *heightArr3,int length3,
							float *corrArr5,float *dbArr5,float *heightArr5,int refLength,
							float light,float *outFre,int *valid,
							int *formatFlag,
							float *fre1,float *fre2,float *fre3,
							float *db1,float *db2,float *db3){

	int flag=0;
	float fre=0;

	float *corrArr2=NULL;
	float *dbArr2=NULL;
	float *heightArr2=NULL;

	int *indexArr1=NULL;
	int *indexArr2=NULL;

	if(!length1){
		return 0;
	}

	corrArr2=__vnew(length1, NULL);
	dbArr2=__vnew(length1, NULL);
	heightArr2=__vnew(length1, NULL);

	__varangei(0, length1, 1, &indexArr1);
	__varangei(0, length1, 1, &indexArr2);

	memcpy(corrArr2, corrArr1, sizeof(float )*length1);
	memcpy(dbArr2, dbArr1, sizeof(float )*length1);
	memcpy(heightArr2, heightArr1, sizeof(float )*length1);

	// dB desc ->corrArr2
	__vcorrsort1(dbArr2, corrArr2 ,NULL,NULL,length1, 1);
	// fre asc ->indexArr1
	__vcorrsort1(corrArr2, NULL ,NULL,indexArr1,length1, 0);

	// fre desc
	memcpy(dbArr2, dbArr1, sizeof(float )*length1);
	__vcorrsort1(dbArr2, corrArr2 ,heightArr2,NULL,length1, 1);

	// fast
	fre=__queue_standard(corrArr1,dbArr1,heightArr1,indexArr1,length1, 
						corrArr3,dbArr3,heightArr3,length3,
						corrArr5,dbArr5,heightArr5,refLength,
						light,valid,
						formatFlag,
						fre1,fre2,fre3,
						db1,db2,db3);

	free(corrArr2);
	free(dbArr2);
	free(heightArr2);

	free(indexArr1);
	free(indexArr2);

	*outFre=fre;
	if(fre){
		flag=6;
	}

	return flag;
}

static int __trist_cut(float *corrArr1,float *dbArr1,float *heightArr1,int length1,
					float *corrArr3,float *dbArr3,float *heightArr3,int length3,
					float *corrArr5,float *dbArr5,float *heightArr5,int refLength,
					float light,float *outFre,int *valid,
					int *formatFlag,
					float *fre1,float *fre2,float *fre3,
					float *db1,float *db2,float *db3){

	int flag=0;
	float fre=0;

	float *corrArr2=NULL;
	float *dbArr2=NULL;
	float *heightArr2=NULL;

	int *indexArr1=NULL;
	int *indexArr2=NULL;

	if(!length1){
		return 0;
	}

	corrArr2=__vnew(length1, NULL);
	dbArr2=__vnew(length1, NULL);
	heightArr2=__vnew(length1, NULL);

	__varangei(0, length1, 1, &indexArr1);
	__varangei(0, length1, 1, &indexArr2);

	memcpy(corrArr2, corrArr1, sizeof(float )*length1);
	memcpy(dbArr2, dbArr1, sizeof(float )*length1);
	memcpy(heightArr2, heightArr1, sizeof(float )*length1);

	// dB desc ->corrArr2
	__vcorrsort1(dbArr2, corrArr2 ,NULL,NULL,length1, 1);
	// fre asc ->indexArr1
	__vcorrsort1(corrArr2, NULL ,NULL,indexArr1,length1, 0);

	// fre desc
	memcpy(dbArr2, dbArr1, sizeof(float )*length1);
	__vcorrsort1(dbArr2, corrArr2 ,heightArr2,NULL,length1, 1);

	// fast
	fre=__queue_cut(corrArr1,dbArr1,heightArr1,indexArr1,length1, 
					corrArr3,dbArr3,heightArr3,length3,
					corrArr5,dbArr5,heightArr5,refLength, 
					light,valid,
					formatFlag,
					fre1,fre2,fre3,
					db1,db2,db3);

	free(corrArr2);
	free(dbArr2);
	free(heightArr2);

	free(indexArr1);
	free(indexArr2);

	*outFre=fre;
	if(fre){
		flag=1;
	}

	return flag;
}

int __trist_fast(float *corrArr1,float *dbArr1,float *heightArr1,int length1,
				float *corrArr3,float *dbArr3,float *heightArr3,int refLength,
				float light,float *outFre,int *valid,
				int *formatFlag,
				float *fre1,float *fre2,float *fre3,
				float *db1,float *db2,float *db3){

	int flag=0;
	float fre=0;

	float *corrArr2=NULL;
	float *dbArr2=NULL;
	float *heightArr2=NULL;

	int *indexArr1=NULL;
	int *indexArr2=NULL;

	if(!length1){
		return 0;
	}

	corrArr2=__vnew(length1, NULL);
	dbArr2=__vnew(length1, NULL);
	heightArr2=__vnew(length1, NULL);

	__varangei(0, length1, 1, &indexArr1);
	__varangei(0, length1, 1, &indexArr2);

	memcpy(corrArr2, corrArr1, sizeof(float )*length1);
	memcpy(dbArr2, dbArr1, sizeof(float )*length1);
	memcpy(heightArr2, heightArr1, sizeof(float )*length1);

	// dB desc ->corrArr2
	__vcorrsort1(dbArr2, corrArr2 ,NULL,NULL,length1, 1);
	// fre asc ->indexArr1
	__vcorrsort1(corrArr2, NULL ,NULL,indexArr1,length1, 0);

	// fre desc
	memcpy(dbArr2, dbArr1, sizeof(float )*length1);
	__vcorrsort1(dbArr2, corrArr2 ,heightArr2,NULL,length1, 1);

	// fast
	fre=__queue_fast(corrArr1,dbArr1,heightArr1,indexArr1,length1,
					corrArr3,dbArr3,heightArr3,refLength,
					light,valid,
					formatFlag,
					fre1,fre2,fre3,
					db1,db2,db3);

	free(corrArr2);
	free(dbArr2);
	free(heightArr2);

	free(indexArr1);
	free(indexArr2);

	*outFre=fre;
	if(fre){
		flag=2;
	}

	return flag;
}

static int __trist(float *corrArr1,float *dbArr1,float *heightArr1,int length1,float light,float *outFre,int *valid,int *status){
	int flag=0;
	float fre=0;

	float *corrArr2=NULL;
	float *dbArr2=NULL;
	float *heightArr2=NULL;

	int *indexArr1=NULL;
	int *indexArr2=NULL;

	if(!length1){
		return 0;
	}

	corrArr2=__vnew(length1, NULL);
	dbArr2=__vnew(length1, NULL);
	heightArr2=__vnew(length1, NULL);

	__varangei(0, length1, 1, &indexArr1);
	__varangei(0, length1, 1, &indexArr2);

	memcpy(corrArr2, corrArr1, sizeof(float )*length1);
	memcpy(dbArr2, dbArr1, sizeof(float )*length1);
	memcpy(heightArr2, heightArr1, sizeof(float )*length1);

	// dB desc ->corrArr2
	__vcorrsort1(dbArr2, corrArr2 ,NULL,NULL,length1, 1);
	// fre asc ->indexArr1
	__vcorrsort1(corrArr2, NULL ,NULL,indexArr1,length1, 0);

	// fre desc
	memcpy(dbArr2, dbArr1, sizeof(float )*length1);
	__vcorrsort1(dbArr2, corrArr2 ,heightArr2,NULL,length1, 1);

	// direct
	fre=__queue_direct(corrArr1,dbArr1,heightArr1,indexArr1,length1,light,valid);
	if(fre){
		flag=3;
	}
	
	// slide
	if(!fre){
		fre=__queue_slide(corrArr1,dbArr1,heightArr1,indexArr1,length1,light,valid,status);
		if(fre){
			flag=4;
		}
	}

	// weak
	if(!fre){
		fre=__queue_weak(corrArr1,dbArr1,heightArr1,indexArr1,length1,light,valid,status);
		if(fre){
			flag=5;
		}
	}

	free(corrArr2);
	free(dbArr2);
	free(heightArr2);

	free(indexArr1);
	free(indexArr2);

	*outFre=fre;
	return flag;
}

static int __arr_maxIndex(float *arr,int length){
	int index=0;

	float value=0;

	value=arr[0];
	for(int i=1;i<length;i++){
		if(value<arr[i]){
			value=arr[i];
			index=i;
		}
	}

	return index;
}









