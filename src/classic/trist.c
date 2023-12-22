// 

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "../util/flux_util.h"

#include "trist.h"

// 0 desc 1 asc
static void __arr_relateSort(float *arr1,float *arr2,float *arr3,int *arr4,int length,int type);
static int __arr_getIndex(int *arr,int length,int value);
static float __sub_valid(float *freArr,float *heightArr,float *dbArr,int length);

static int __isEqual(float v1,float v2);

int trist(float *corrArr,float *dbArr,float *heightArr,int *midiArr1,int length,
		float *freArr,float *dbArr2,float *heightArr2,int *midiArr2,int count1,int count2,
		float *outFre){
	int flag=0;

	float fre1=0;
	float fre2=0;
	float fre3=0;

	float sub1=0;
	float sub2=0;
	float sub3=0;

	int k1=0;
	int k2=0;
	int k3=0;
	int k4=0;
	int k5=0;

	float *arr1=NULL;
	float *arr2=NULL;

	float _fre=0;

	sub1=fabsf(corrArr[0]-corrArr[1]);
	sub2=fabsf(corrArr[0]-corrArr[2]);
	sub3=fabsf(corrArr[1]-corrArr[2]);
	
	arr1=__vnew(12, NULL);
	arr2=__vnew(12, NULL);

	// 123!!!
	if(!flag){
		__vsort(corrArr, 3, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		if(k1==2&&k2==3){

			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 1234!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		if(k1==2&&k2==3&&k3==4){
			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 1234nn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		k4=util_freTimes(arr1[4], arr1[0]);
		k5=util_freTimes(arr1[5], arr1[0]);
		if(k1==2&&k2==3&&k3==4&&
			k4&&k5){
			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 1234n!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		k4=util_freTimes(arr1[4], arr1[0]);
		if(k1==2&&k2==3&&k3==4&&
			k4){
			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 1247!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		if(k1==2&&k2==4&&k3==7){
			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 1234*!!!
	{	
		// 1*234!!!
		if(!flag){
			__vsort(corrArr, 5, 0, arr1);
			k1=util_freTimes(arr1[2], arr1[0]);
			k2=util_freTimes(arr1[3], arr1[0]);
			k3=util_freTimes(arr1[4], arr1[0]);
			if((k1==2&&k2==3&&k3==4)&&
				!__isEqual(arr1[1], corrArr[0])){

				_fre=arr1[2]/2;
				flag=1;
			}
		}

		// 12*34!!!
		if(!flag){
			__vsort(corrArr, 5, 0, arr1);
			k1=util_freTimes(arr1[1], arr1[0]);
			k2=util_freTimes(arr1[3], arr1[0]);
			k3=util_freTimes(arr1[4], arr1[0]);
			if((k1==2&&k2==3&&k3==4)&&
				!__isEqual(arr1[2], corrArr[0])){

				_fre=arr1[1]/2;
				flag=1;
			}
		}

		// 123*4!!!
		if(!flag){
			__vsort(corrArr, 5, 0, arr1);
			k1=util_freTimes(arr1[1], arr1[0]);
			k2=util_freTimes(arr1[2], arr1[0]);
			k3=util_freTimes(arr1[4], arr1[0]);
			if((k1==2&&k2==3&&k3==4)&&
				!__isEqual(arr1[3], corrArr[0])){

				_fre=arr1[1]/2;
				flag=1;
			}
		}

		// 1234*!!!
		if(!flag){
			__vsort(corrArr, 5, 0, arr1);
			k1=util_freTimes(arr1[1], arr1[0]);
			k2=util_freTimes(arr1[2], arr1[0]);
			k3=util_freTimes(arr1[3], arr1[0]);
			if((k1==2&&k2==3&&k3==4)&&
				!__isEqual(arr1[4], corrArr[0])){

				_fre=arr1[1]/2;
				flag=1;
			}
		}

		// *1234!!!
		if(!flag){
			__vsort(corrArr, 5, 0, arr1);
			k1=util_freTimes(arr1[2], arr1[1]);
			k2=util_freTimes(arr1[3], arr1[1]);
			k3=util_freTimes(arr1[4], arr1[1]);
			if((k1==2&&k2==3&&k3==4)&&
				!__isEqual(arr1[0], corrArr[0])){

				_fre=arr1[2]/2;
				flag=1;
			}
		}
	}

	// 123*!!!
	{
		// 1*23
		if(!flag){
			__vsort(corrArr, 4, 0, arr1);
			k1=util_freTimes(arr1[2], arr1[0]);
			k2=util_freTimes(arr1[3], arr1[0]);
			if((k1==2&&k2==3)&&
				!__isEqual(arr1[1], corrArr[0])){

				_fre=arr1[2]/2;
				flag=1;
			}
		}

		// 12*3
		if(!flag){
			__vsort(corrArr, 4, 0, arr1);
			k1=util_freTimes(arr1[1], arr1[0]);
			k2=util_freTimes(arr1[3], arr1[0]);
			if((k1==2&&k2==3)&&
				!__isEqual(arr1[2], corrArr[0])){

				_fre=arr1[1]/2;
				flag=1;
			}
		}

		// 123*
		if(!flag){
			__vsort(corrArr, 4, 0, arr1);
			k1=util_freTimes(arr1[1], arr1[0]);
			k2=util_freTimes(arr1[2], arr1[0]);
			if((k1==2&&k2==3)&&
				!__isEqual(arr1[3], corrArr[0])){

				_fre=arr1[1]/2;
				flag=1;
			}
		}

		// *123
		if(!flag){
			__vsort(corrArr, 4, 0, arr1);
			k1=util_freTimes(arr1[2], arr1[1]);
			k2=util_freTimes(arr1[3], arr1[1]);
			if((k1==2&&k2==3)&&
				!__isEqual(arr1[0], corrArr[0])){

				_fre=arr1[2]/2;
				flag=1;
			}
		}
	}

	// 1*23nn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[2], arr1[0]);
		k2=util_freTimes(arr1[3], arr1[0]);
		k3=util_freTimes(arr1[4], arr1[0]);
		k4=util_freTimes(arr1[5], arr1[0]);
		if(k1==2&&k2==3&&k3&&k4&&
			!__isEqual(arr1[1], corrArr[0])){

			_fre=arr1[2]/2;
			flag=1;
		}
	}

	// 123nn!!!
	if(!flag){
		__vsort(corrArr, 5, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		k4=util_freTimes(arr1[4], arr1[0]);
		if(k1==2&&k2==3&&k3&&k4){
			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 123nnn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		k4=util_freTimes(arr1[4], arr1[0]);
		k5=util_freTimes(arr1[5], arr1[0]);
		if(k1==2&&k2==3&&k3&&k4&&k5){
			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 123n!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		if(k1==2&&k2==3&&k3){
			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 12!!! +dB
	// if(!flag){
	// 	__vsort(corrArr, 2, 0, arr1);
	// 	k1=util_freTimes(arr1[1], arr1[0]);
	// 	if(k1==2&&
	// 		roundf(dbArr[0]-dbArr[2])>=10){

	// 		_fre=arr1[1]/2;
	// 		flag=1;
	// 	}
	// }

	// // 12 +dB
	// if(!flag){
	// 	k1=util_freTimes(corrArr[1], corrArr[0]);
	// 	if(k1==2&&
	// 		corrArr[1]>corrArr[0]&&
	// 		roundf(dbArr[0]-dbArr[2])>=10){

	// 		_fre=corrArr[0];
	// 		flag=1;
	// 	}
	// }

	// // 13 +dB
	// if(!flag){
	// 	k1=util_freTimes(corrArr[1], corrArr[0]);
	// 	if(k1==3&&
	// 		corrArr[1]>corrArr[0]&&
	// 		roundf(dbArr[0]-dbArr[1])>=10){

	// 		_fre=corrArr[0];
	// 		flag=1;
	// 	}
	// }

	// 1? +23!!! +dB
	if(!flag){
		__vsort(corrArr+2, 2, 0, arr1);
		k1=util_freTimes(arr1[0], corrArr[0]);
		k2=util_freTimes(arr1[1], corrArr[0]);
		// k3=util_freTimes(corrArr[1], corrArr[0]);
		if(k1==2&&k2==3&&
			arr1[0]>corrArr[0]&&
			roundf(dbArr[0]-dbArr[1])>=10&&
			roundf(dbArr[0]-dbArr[2])>=10&&
			roundf(dbArr[0]-dbArr[3])>=10){

			_fre=corrArr[0];
			flag=1;
		}
	}

	// 1 +dB
	if(!flag){
		if(roundf(fabs(dbArr[0]))>=48&&
			roundf(dbArr[0]-dbArr[1])>=20){

			_fre=corrArr[0];
			flag=1;
		}
	}

	// 12468!!!-->246/468 height连续
	if(!flag){
		__vsort(corrArr, 5, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		k4=util_freTimes(arr1[4], arr1[0]);
		if(k1==2&&k2==4&&k3==6&&k4==8){
			_fre=arr1[2]/2;
			flag=1;
		}
	}

	// 1246!!!+dB
	if(!flag){
		memcpy(arr1, corrArr, sizeof(float )*4);
		memcpy(arr2, dbArr, sizeof(float )*4);
		__arr_relateSort(arr1,arr2,NULL,NULL,4,1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		if((k1==2&&k2==4&&k3==6)){
			if(dbArr[0]-arr2[0]<6){
				_fre=arr1[1]/2;
				
			}
			else{
				_fre=arr1[1];
			}
			flag=1;
		}
	}

	// // 12367!!!/1367!!!
	// if(!flag){
	// 	__vsort(corrArr, 5, 0, arr1);
	// 	k1=util_freTimes(arr1[1], arr1[0]);
	// 	k2=util_freTimes(arr1[2], arr1[0]);
	// 	k3=util_freTimes(arr1[3], arr1[0]);
	// 	k4=util_freTimes(arr1[4], arr1[0]);
	// 	if(k1==2&&k2==3&&k3==6&&k4==7){
	// 		_fre=arr1[1]/2;
	// 		flag=1;
	// 	}

	// 	if(!flag){
	// 		__vsort(corrArr, 4, 0, arr1);
	// 		k1=util_freTimes(arr1[1], arr1[0]);
	// 		k2=util_freTimes(arr1[2], arr1[0]);
	// 		k3=util_freTimes(arr1[3], arr1[0]);
	// 		if(k1==3&&k2==6&&k3==7){
	// 			_fre=arr1[1]/3;
	// 			flag=1;
	// 		}
	// 	}
	// }

	// 125nn!!!
	if(!flag){
		__vsort(corrArr, 5, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		k4=util_freTimes(arr1[4], arr1[0]);
		if(k1==2&&k2==5&&
			k3&&k4){

			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 12457!!!
	if(!flag){
		__vsort(corrArr, 5, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		k4=util_freTimes(arr1[4], arr1[0]);
		if(k1==2&&k2==4&&k3==5&&k4==7){

			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 2357!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		if(k1==3&&k2==5&&k3==7){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 2367!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		if(k1==3&&k2==6&&k3==7){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 2347!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		if(k1==3&&k2==4&&k3==7){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 2346!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		if(k1==3&&k2==4&&k3==6){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 237n!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		if(k1==3&&k2==7&&k3){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 237nn!!!
	if(!flag){
		__vsort(corrArr, 5, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		if(k1==3&&k2==7&&k3&&k4){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 2367n!!!
	if(!flag){
		__vsort(corrArr, 5, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		if(k1==3&&k2==6&&k3==7&&k4){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 124[5|7|9|11]nn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]);
		k2=util_freTimes(arr1[2], arr1[0]);
		k3=util_freTimes(arr1[3], arr1[0]);
		k4=util_freTimes(arr1[4], arr1[0]);
		k5=util_freTimes(arr1[5], arr1[0]);
		if(k1==2&&k2==4&&
			k3&&k4&&k5&&
			((k3==5||k3==7||k3==9||k3==11)||
			(k4==5||k4==7||k4==9||k4==11)||
			(k5==5||k5==7||k5==9||k5==11)) ){

			_fre=arr1[1]/2;
			flag=1;
		}
	}

	// 234nnn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		k5=util_freTimes(arr1[5], arr1[0]/2);
		if(k1==3&&k2==4&&
			k3&&k4&&k5){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 234nn!!!
	if(!flag){
		__vsort(corrArr, 5, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		if(k1==3&&k2==4&&
			k3&&k4){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 234*
	{
		// 2*34
		if(!flag){
			__vsort(corrArr, 4, 0, arr1);
			k1=util_freTimes(arr1[2], arr1[0]/2);
			k2=util_freTimes(arr1[3], arr1[0]/2);
			if((k1==3&&k2==4)&&
				!__isEqual(arr1[1], corrArr[0])){

				_fre=arr1[0]/2;
				flag=1;
			}
		}

		// 23*4
		if(!flag){
			__vsort(corrArr, 4, 0, arr1);
			k1=util_freTimes(arr1[1], arr1[0]/2);
			k2=util_freTimes(arr1[3], arr1[0]/2);
			if((k1==3&&k2==4)&&
				!__isEqual(arr1[2], corrArr[0])){

				_fre=arr1[0]/2;
				flag=1;
			}
		}

		// 234*
		if(!flag){
			__vsort(corrArr, 4, 0, arr1);
			k1=util_freTimes(arr1[1], arr1[0]/2);
			k2=util_freTimes(arr1[2], arr1[0]/2);
			if((k1==3&&k2==4)&&
				!__isEqual(arr1[3], corrArr[0])){

				_fre=arr1[0]/2;
				flag=1;
			}
		}

		// *234
		if(!flag){
			__vsort(corrArr, 4, 0, arr1);
			k1=util_freTimes(arr1[2], arr1[1]/2);
			k2=util_freTimes(arr1[3], arr1[1]/2);
			if((k1==3&&k2==4)&&
				!__isEqual(arr1[0], corrArr[0])){

				_fre=arr1[1]/2;
				flag=1;
			}
		}
	}

	// 245nnn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		k5=util_freTimes(arr1[5], arr1[0]/2);
		if(k1==4&&k2==5&&
			k3&&k4&&k5){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 245n!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		if(k1==4&&k2==5&&k3){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 234n!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		if(k1==3&&k2==4){
			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 257nnn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		k5=util_freTimes(arr1[5], arr1[0]/2);
		if(k1==5&&k2==7&&
			k3&&k4&&k5){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 267nnn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		k5=util_freTimes(arr1[5], arr1[0]/2);
		if(k1==6&&k2==7&&
			k3&&k4&&k5){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 279nnn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		k5=util_freTimes(arr1[5], arr1[0]/2);
		if(k1==7&&k2==9&&
			k3&&k4&&k5){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 2467nn!!!
	if(!flag){
		__vsort(corrArr, 6, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		k5=util_freTimes(arr1[5], arr1[0]/2);
		if(k1==4&&k2==6&&k3==7&&
			k4&&k5){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 234nn!!!
	if(!flag){
		__vsort(corrArr, 5, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		k4=util_freTimes(arr1[4], arr1[0]/2);
		if(k1==3&&k2==4&&
			k3&&k4){

			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 234n!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/2);
		k2=util_freTimes(arr1[2], arr1[0]/2);
		k3=util_freTimes(arr1[3], arr1[0]/2);
		if(k1==3&&k2==4&&k3){
			_fre=arr1[0]/2;
			flag=1;
		}
	}

	// 3456!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/3);
		k2=util_freTimes(arr1[2], arr1[0]/3);
		k3=util_freTimes(arr1[3], arr1[0]/3);
		if(k1==4&&k2==5&&k3==6){
			_fre=arr1[0]/3;
			flag=1;
		}
	}

	// 3467!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/3);
		k2=util_freTimes(arr1[2], arr1[0]/3);
		k3=util_freTimes(arr1[3], arr1[0]/3);
		if(k1==4&&k2==6&&k3==7){
			_fre=arr1[0]/3;
			flag=1;
		}
	}

	// 3567!!!
	if(!flag){
		__vsort(corrArr, 4, 0, arr1);
		k1=util_freTimes(arr1[1], arr1[0]/3);
		k2=util_freTimes(arr1[2], arr1[0]/3);
		k3=util_freTimes(arr1[3], arr1[0]/3);
		if(k1==5&&k2==6&&k3==7){
			_fre=arr1[0]/3;
			flag=1;
		}
	}

	// 32n
	if(!flag){
		if(corrArr[0]>corrArr[1]&&
			fabsf(dbArr[1]-dbArr[2])<6){
			k1=util_freTimes(sub1, corrArr[0]);
			k2=util_freTimes(sub1, corrArr[1]);
			if(k1==3&&k2==2){
				_fre=corrArr[1]/2;
				flag=1;
			}
		}
	}

	// nnn!!!
	if(!flag){
		int midi=0;
		int index=0;

		midi=util_freToMidi(corrArr[0]);
		index=__arr_getIndex(midiArr2, count1+count2, midi);
		if(index!=-1){
			if(index-1>=0){
				fre1=freArr[index-1];
				sub1=fabsf(fre1-corrArr[0]);
				k1=util_freTimes(sub1, corrArr[0]);
				k2=util_freTimes(sub1, corrArr[1]);
				k3=util_freTimes(sub1, corrArr[2]);
				if(k1&&k2&&k3){
					_fre=corrArr[0]/k1;
					flag=1;
				}
			}

			if(!flag&&index+1<count1+count2){
				fre1=freArr[index+1];
				sub1=fabsf(fre1-corrArr[0]);
				k1=util_freTimes(sub1, corrArr[0]);
				k2=util_freTimes(sub1, corrArr[1]);
				k3=util_freTimes(sub1, corrArr[1]);
				if(k1&&k2&&k3){
					_fre=corrArr[0]/k1;
					flag=1;
				}
			}
		}
	}

	// 1nn
	if(!flag){
		fre1=corrArr[0];
		if(corrArr[1]>corrArr[0]&&
			corrArr[2]>corrArr[0]){
			k1=util_freTimes(fre1, corrArr[1]);
			k2=util_freTimes(fre1, corrArr[2]);
			if(k1&&k2){
				_fre=corrArr[1]/k1;
				flag=1;
			}
		}
	}

	// 2nn
	if(!flag){
		fre1=corrArr[0]/2;
		if(corrArr[1]>corrArr[0]&&
			corrArr[2]>corrArr[0]){
			k1=util_freTimes(fre1, corrArr[1]);
			k2=util_freTimes(fre1, corrArr[2]);
			if(k1&&k2){
				_fre=fre1;
				flag=1;
			}
		}
	}

	// n2n
	if(!flag){
		fre1=corrArr[1]/2;
		if(corrArr[0]>corrArr[1]&&
			corrArr[2]>corrArr[1]){
			k1=util_freTimes(fre1, corrArr[0]);
			k2=util_freTimes(fre1, corrArr[2]);
			if(k1&&k2){
				_fre=fre1;
				flag=1;
			}
		}
	}

	// 23
	if(!flag){ 
		k1=util_freTimes(corrArr[1], corrArr[0]/2);
		if(k1==3&&
			corrArr[0]<corrArr[1]){

			_fre=corrArr[0]/2;
			flag=1;
		}
	}

	*outFre=_fre;

	free(arr1);
	free(arr2);
	return flag;
}

// 0 desc 1 asc
static void __arr_relateSort(float *arr1,float *arr2,float *arr3,int *arr4,int length,int type){
	float value1=0;
	float value2=0;
	float value3=0;
	int value4=0;

	for(int i=0;i<length;i++){
		for(int j=i+1;j<length;j++){
			if(type==0){ // 降序
				if(arr1[i]<arr1[j]){
					value1=arr1[i];
					arr1[i]=arr1[j];
					arr1[j]=value1;

					if(arr2){
						value2=arr2[i];
						arr2[i]=arr2[j];
						arr2[j]=value2;
					}
					
					if(arr3){
						value3=arr3[i];
						arr3[i]=arr3[j];
						arr3[j]=value3;
					}

					if(arr4){
						value4=arr4[i];
						arr4[i]=arr4[j];
						arr4[j]=value4;
					}
				}
			}
			else{ // 升序
				if(arr1[i]>arr1[j]){
					value1=arr1[i];
					arr1[i]=arr1[j];
					arr1[j]=value1;

					if(arr2){
						value2=arr2[i];
						arr2[i]=arr2[j];
						arr2[j]=value2;
					}
					
					if(arr3){
						value3=arr3[i];
						arr3[i]=arr3[j];
						arr3[j]=value3;
					}

					if(arr4){
						value4=arr4[i];
						arr4[i]=arr4[j];
						arr4[j]=value4;
					}
				}
			}
		}
	}
}

static int __arr_getIndex(int *arr,int length,int value){
	int index=-1;

	for(int i=0;i<length;i++){
		if(arr[i]==value){
			index=i;
			break;
		}
	}

	return index;
}

static float __sub_valid(float *freArr,float *heightArr,float *dbArr,int length){
	float value=0;



	return value;
}

static int __isEqual(float v1,float v2){
	int flag=0;

	if(fabsf(v1-v2)<0.1){
		flag=1;
	}

	return flag;
}

