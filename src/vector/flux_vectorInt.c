// 

#include <string.h>
#include <math.h>

#include "flux_vector.h"

int *__vnewi(int length,int *value){
	int *arr=NULL;
	int _value=0;

	arr=(int *)calloc(length, sizeof(int ));
	if(value){
		_value=*value;
	}

	if(_value!=0){
		for(int i=0;i<length;i++){
			arr[i]=_value;
		}
	}
	
	return arr;
}

int __varangei(int start,int stop,int step,int **outArr){
	int *arr=NULL;
	int length=0;

	if(stop<=start||!outArr){
		return 0;
	}

	length=ceilf((stop-start)/step);
	arr=__vnewi(length, NULL);
	for(int i=0;i<length;i++){
		arr[i]=start+i*step;
	}

	*outArr=arr;
	return length;
}

void __vfilli(int *arr,int length,int value){

	for(int i=0;i<length;i++){
		arr[i]=value;
	}
}

void __vcopyi(int *dstArr,int *srcArr,int length){

	memcpy(dstArr, srcArr, sizeof(float )*length);
}

void __vdebugi(int *vArr1,int length,int type){
	if(type){
		printf("vector is:\n");

		printf("	");
		for(int i=0;i<length;i++){
			printf("%d:%d ,",i,vArr1[i]);
		}
	}
	else{
		printf("vector([");
		for(int i=0;i<length-1;i++){
			printf("%d, ",vArr1[i]);
		}
		printf("%d",vArr1[length-1]);
		printf("])\n");
	}
}

void __mdebugi(int *mArr1,int nLength,int mLength,int type){
	if(type){
		printf("matrix is:\n");

		for(int i=0;i<nLength;i++){
			printf("	%d:\n",i);
			printf("		");
			for(int j=0;j<mLength;j++){
				printf("%d:%d ,",j,mArr1[i*mLength+j]);
			}

			printf("\n\n");
		}
	}
	else{
		printf("matrix([");
		for(int i=0;i<nLength-1;i++){
			printf("[");

			for(int j=0;j<mLength-1;j++){
				printf("%d ,",mArr1[i*mLength+j]);
			}
			printf("%d",mArr1[i*mLength+mLength-1]);

			printf("],");
			printf("\n        ");
		}

		printf("[");
		for(int j=0;j<mLength-1;j++){
			printf("%d ,",mArr1[(nLength-1)*mLength+j]);
		}
		printf("%d",mArr1[(nLength-1)*mLength+mLength-1]);

		printf("]");
		printf("])\n");
	}
}

void __vsubi(int *vArr1,int *vArr2,int length,int *vArr3){
	int *arr=NULL;

	if(vArr3){
		arr=vArr3;
	}
	else{
		arr=vArr1;
	}

	for(int i=0;i<length;i++){
		arr[i]=vArr1[i]-vArr2[i];
	}
}


