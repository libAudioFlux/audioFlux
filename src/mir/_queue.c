// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "_queue.h"

typedef struct{
	int *indexArr;
	int *kArr;

	int *numArr;

	int length;
	int capLength;

} QueueMap;

static float __queue_weakValid(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length);
static float __queue_cutValid(float *freArr,float *dbArr,int length,int flag,int mode,
							float *freArr2,float *dbArr2,int length2,
							float *freArr3,float *dbArr3,int length3);

static int __queue_query(float *freArr,float *dbArr,float *heightArr,int length,float value);

static int __queue_query2(float *freArr,float *dbArr,float *heightArr,int length,int start,float value1,float value2,int strict);
static int __queue_query3(float *freArr,float *dbArr,float *heightArr,int length,int start,float value1,float value2,int strict);

static int __queue_query2Inf(float *freArr,float *dbArr,float *heightArr,int length,int start,float value1,float value2,int strict);
static int __queue_four(float *freArr,float *dbArr,float *heightArr,int length,float value);
static int __queue_valid98(float *freArr,float *dbArr,float *heightArr,int length,int start,float value,int strict);
static int __queue_odd98(float *freArr,float *dbArr,float *heightArr,int length,int start,float value);

static float __queue_slideValid(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,float value);

// one/two/jump, not move
static int __queue_cal(float *freArr,float *dbArr,float *heightArr,int length,int start,int flag,
					int *index1,int *k1,
					int *index2,int *k2,
					int *jumpFlag);

static float __queue_twoMove(float *freArr,float *dbArr,float *heightArr,int length,int start,
							int index1,int k1,
							int index2,int k2,
							int jumpFlag,
							int *offset);

static float __queue_oneMove(float *freArr,float *dbArr,float *heightArr,int length,int start,
							int index1,int k1,
							int *index2,int *k2,
							int *offset);

static float __queue_jumpMove(float *freArr,float *dbArr,float *heightArr,int length,int start,
							int index1,int k1,
							int jumpFlag,
							int *index2,int *k2,
							int *offset);

static float __queue_jumpBound(float *freArr,float *dbArr,float *heightArr,int length,int start,
							int index1,int k1,
							int jumpFlag,
							int *index2,int *k2,
							int *offset);

static int __queue_has(float *freArr,int length,float baseFre,int start,int *index);
static int __queue_isEqual(float *freArr,int length,int index1,int k1,int index2,int k2);

static float __checkFre(float fre1,float fre2,float fre3);
static float __checkFre2(float fre1,float fre2,float fre3);

static int __arr_getIndex(int *arr,int length,int value);
static int __arr_maxIndex(float *arr,int length);
static int __arr_cut(float *arr,int length,float value);

static int __validFre3(float fre1,float fre2,float fre3,float base,int k1,int k2,int k3);
static int __validFre2(float fre1,float fre2,float base,int k1,int k2);
static int __isValidTimes(int k1,int k2,int k3);

static void __map_add(QueueMap *map,float *freArr,int length,int index,int k);
static int __map_has(QueueMap *map,float *freArr,int length,int index,int k);

static QueueMap *__createMap(int capLength);
static void __freeMap(QueueMap *queueMap);

/***
	1:1/1:2/1:3 2:2 2:3
	1:4 ???
****/
float __queue_fre3(float value1,float value2,float value3,
					int *s1,int *s2,
					int *k1,int *k2,int *k3){
	float base=0;

	float sub1=0;
	float sub2=0;
	float _sub=0;

	int _k=0;
	int _k1=0;
	int _k2=0;
	int _k3=0;

	int _s1=0;
	int _s2=0;

	int gFlag=0;
	int type=0;

	sub1=value2-value1;
	sub2=value3-value2;
	if(sub1>sub2){
		_sub=sub1;
		sub1=sub2;
		sub2=_sub;

		gFlag=1;
	}

	_k=util_calRangeTimes(sub1,sub2,&type);
	if(_k==1){ // 1:1
		_k1=util_calRangeTimes(sub1,value1,&type);
		_k2=util_calRangeTimes(sub1,value2,&type);
		if(_k1&&_k2){ 
			// _k3=util_calRangeTimes(sub1,value3,&type);
			_k3=_k2+1;

			_s1=1;
			_s2=1;

			base=value1/_k1;
		}
		else{ // 2:2
			_k1=util_calRangeTimes(sub1/2,value1,&type);
			_k2=util_calRangeTimes(sub1/2,value2,&type);
			if(_k1&&_k2){
				// _k3=util_calRangeTimes(sub1/2,value3,&type);
				_k3=_k2+2;

				if(_k1%2==1){ // ->2:2
					_s1=2;
					_s2=2;

					base=value1/_k1;
				}
				else{ // ->1:1
					_s1=1;
					_s2=1;

					_k1/=2;
					_k2/=2;
					_k3/=2;

					base=value1/_k1;
				}
			}
		}
	}
	else if(_k>=2&&_k<=4){ // 1:2 1:3 1:4
		_k1=util_calRangeTimes(sub1,value1,&type);
		_k2=util_calRangeTimes(sub1,value2,&type);
		if(_k1&&_k2){
			// _k3=util_calRangeTimes(sub1,value3,&type);
			if(gFlag){
				_k3=_k2+1;
			}
			else{
				_k3=_k2+_k;
			}

			_s1=(gFlag?_k:1);
			_s2=(gFlag?1:_k);

			base=value1/_k1;
		}
	}
	else{ // 2:3
		_sub=sub2-sub1;
		_k1=util_calRangeTimes(_sub,sub1,&type);
		_k2=util_calRangeTimes(_sub,sub2,&type);
		if(_k1==2&&_k2==3){
			_k1=util_calRangeTimes(sub1/2,value1,&type);
			_k2=util_calRangeTimes(sub1/2,value2,&type);
			if(_k1&&_k2){
				// _k3=util_calRangeTimes(sub1/2,value3,&type);
				if(gFlag){
					_k3=_k2+2;
				}
				else{
					_k3=_k2+3;
				}

				_s1=(gFlag?3:2);
				_s2=(gFlag?2:3);

				base=value1/_k1;
			}
		}
	}

	if(!base){
		_k=roundf(sub2/sub1);
		if(_k==1){ // 1:1
			_k1=roundf(value1/sub1);
			_k2=roundf(value2/sub1);
			if(_k1+1==_k2){ 
				// _k3=util_calRangeTimes(sub1,value3,&type);
				_k3=_k2+1;

				_s1=1;
				_s2=1;

				base=value1/_k1;
			}
			else{ // 2:2
				_k1=roundf(value1/(sub1/2));
				_k2=roundf(value2/(sub1/2));
				if(_k1+2==_k2){
					// _k3=util_calRangeTimes(sub1/2,value3,&type);
					_k3=_k2+2;

					_s1=2;
					_s2=2;

					base=value1/_k1;
				}
			}
		}
		else if(_k>=2&&_k<=4){ // 1:2 1:3 1:4
			_k1=roundf(value1/sub1);
			_k2=roundf(value2/sub1);
			if(_k1&&_k2){
				// _k3=util_calRangeTimes(sub1,value3,&type);
				if(gFlag){
					_k3=_k2+1;
				}
				else{
					_k3=_k2+_k;
				}

				_s1=(gFlag?_k:1);
				_s2=(gFlag?1:_k);

				base=value1/_k1;
			}
		}

		if(base){
			if(!(fabsf(value2-value1/_k1*_k2)<5&&
				fabsf(value3-value1/_k1*_k3)<5)){

				base=0;
			}
		}
	}

	if(base){
		if(s1){
			*s1=_s1;
		}
		if(s2){
			*s2=_s2;
		}

		if(k1){
			*k1=_k1;
		}
		if(k2){
			*k2=_k2;
		}
		if(k3){
			*k3=_k3;
		}
	}
	else{
		if(s1){
			*s1=0;
		}
		if(s2){
			*s2=0;
		}

		if(k1){
			*k1=0;
		}
		if(k2){
			*k2=0;
		}
		if(k3){
			*k3=0;
		}
	}

	return base;
}

/***
	k=value2/value1
	1:n/2:n
****/
float __queue_fre2(float value1,float value2,
					int *k1,int *k2){
	float fre=0;

	int type=0;
	
	int k=0;
	int _k1=0;
	int _k2=0;

	float sub=0;

	int flag=0;

	// k=util_calFreTimes(value1,value2,NULL);
	k=util_calRangeTimes(value1,value2,NULL);
	if(k){ // >=1
		fre=value1;

		_k1=1;
		_k2=k;
	}
	else{
		sub=value2-value1;

		_k2=util_calRangeTimes(sub,value2,NULL);
		if(_k2){
			_k1=util_calRangeTimes(sub,value1,&type);
			if(_k1&&!type){
				fre=value1/_k1;

				flag=1;
			}
		}

		if(!flag){
			sub/=2;
			_k2=util_calRangeTimes(sub,value2,NULL);
			if(_k2){
				_k1=util_calRangeTimes(sub,value1,&type);
				if(_k1&&!type){
					fre=value1/_k1;
				}
			}
		}
	}

	if(fre){
		if(k1){
			*k1=_k1;
		}
		if(k2){
			*k2=_k2;
		}
	}
	else{
		if(k1){
			*k1=0;
		}
		if(k2){
			*k2=0;
		}
	}

	return fre;
}

int __queue_num(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,int **indexArr2,int **kArr,int **numArr){
	int qLen=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int *_indexArr=NULL;
	int *_kArr=NULL;
	int *_numArr=NULL;

	if(length<3){
		return 0;
	}

	_indexArr=__vnewi(length, NULL);
	_kArr=__vnewi(length, NULL);
	_numArr=__vnewi(length, NULL);
	for(int i=0;i<length-2;i++){
		float _fre=0;

		int _flag=0;
		int _index=0;

		_fre=__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
						&us1,&us2,
						&uk1,&uk2,&uk3);
		if(_fre){
			for(int j=0;j<qLen;j++){
				float _k=0;

				_index=_indexArr[j];
				_k=util_calRangeTimes(freArr[_index]/_kArr[j],_fre,NULL);
				if(_k==1){
					_flag=1;
					_index=j;
					break;
				}
			}

			if(_flag){ // has
				_numArr[_index]++;
			}
			else{
				_indexArr[qLen]=i;
				_kArr[qLen]=uk1;
				_numArr[qLen]=1;

				qLen++;
			}
		}
	}

	if(qLen){
		if(indexArr2){
			*indexArr2=_indexArr;
		}
		else{
			free(_indexArr);
		}

		if(kArr){
			*kArr=_kArr;
		}
		else{
			free(_kArr);
		}

		if(numArr){
			*numArr=_numArr;
		}
		else{
			free(_numArr);
		}
	}
	else{
		free(_indexArr);
		free(_kArr);
		free(_numArr);
	}

	return qLen;
}

/***
	num >=1
	subType 0,1:1+1:2/2:1+2:2; 1,1:1+2:2
	unionType 0, 0; 1, has 1; 2, has 2
	direction 0,start; 1,end
****/
float __queue_multi(float *freArr,float *dbArr,float *heightArr,int length,int num,int subType,int unionType,int direction){
	float fre=0;

	int qLen=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int *_indexArr=NULL;
	int *_kArr=NULL;
	int *_numArr=NULL;

	int step=0;

	if(length<5||num<1){
		return 0;
	}

	if(!unionType){
		step=2;
	}
	else if(unionType==1){
		step=1;
	}

	_indexArr=__vnewi(length, NULL);
	_kArr=__vnewi(length, NULL);
	_numArr=__vnewi(length, NULL);
	if(!direction){
		for(int i=(!direction?0:length-1);(!direction?i<length-2:i>1);i=i+(!direction?1:-1)){
			float _fre=0;

			int _flag=0;
			int _index=0;

			int _sFlag=0;

			_fre=__queue_fre3(freArr[i],freArr[i+(!direction?1:-1)],freArr[i+(!direction?2:-2)],
							&us1,&us2,
							&uk1,&uk2,&uk3);

			if(!subType){ // 1:1+1:2/2:1+2:2
				if((us1==1||us1==2)&&
					(us2==1||us2==2)){

					_sFlag=1;
				}
			}	
			else{
				if(us1==us2&&(us1==1||us1==2)){ // 1:1+2:2
					_sFlag=1;
				}
			}

			if(_fre&&_sFlag){
				for(int j=0;j<qLen;j++){
					float _k=0;

					_index=_indexArr[j];
					_k=util_calRangeTimes(freArr[_index]/_kArr[j],_fre,NULL);
					if(_k==1){
						_flag=1;
						_index=j;
						break;
					}
				}

				if(_flag){ // has
					_numArr[_index]++;

					if(_numArr[_index]==num){
						fre=_fre;
						break;
					}
				}
				else{
					_indexArr[qLen]=i;
					_kArr[qLen]=uk1;
					_numArr[qLen]=1;

					qLen++;

					if(_numArr[qLen]==num){
						fre=_fre;
						break;
					}
				}

				i=i+(!direction?step:-step);
			}
		}
	}
	
	free(_indexArr);
	free(_kArr);
	free(_numArr);

	return fre;
}

int __queue_bear(float *freArr,float *dbArr,float *heightArr,int length,float min,float base,int *index){
	int flag=0;

	float fre=0;
	int k1=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int start=0;

	if(index){
		if(*index>=0){
			start=*index;
		}
	}

	for(int i=start;i<length-2;i++){
		if(freArr[i]>min){
			fre=__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
							&us1,&us2,
							&uk1,&uk2,&uk3);

			if(fre&&
				(us1==1||us1==2)&&
				(us2==1||us2==2)){

				k1=util_calRangeTimes(fre,base,NULL);
				if(k1==1){
					flag=1;
					if(index){
						*index=i;
					}

					break;
				}
			}
		}
	}

	return flag;
}

int __queue_count(float *freArr,float *dbArr,float *heightArr,int length,int start,float min,float base,int step){
	int count=0;

	float fre=0;
	int k1=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	for(int i=start;i<length-2;i++){
		if(freArr[i]>min){
			fre=__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
							&us1,&us2,
							&uk1,&uk2,&uk3);

			if(fre&&
				(us1==1||us1==2)&&
				(us2==1||us2==2)){

				k1=util_calRangeTimes(fre,base,NULL);
				if(k1==1){
					count++;
					i+=step;
				}
			}
		}
	}

	return count;
}

static int __queue_query(float *freArr,float *dbArr,float *heightArr,int length,float value){
	int flag=0;

	for(int i=0;i<length;i++){
		int k=0;
		int type=0;

		k=util_calRangeTimes(value, freArr[i], &type);
		if(k&&!type){
			if(k==5||k==7||k==9||k==11||k==13){
				if(k==5&&fabsf(value*k-freArr[i])<6){
					flag=1;
				}
				else if(k<=9&&fabsf(value*k-freArr[i])<12){
					flag=1;
				}
				else if(k==11&&fabsf(value*k-freArr[i])<18){
					flag=1;
				}
				else if(k==13&&fabsf(value*k-freArr[i])<20){
					flag=1;
				}
			}

			if(flag){
				break;
			}
		}
	}

	return flag;
}

static int __queue_four(float *freArr,float *dbArr,float *heightArr,int length,float value){
	int flag=0;

	int k1=0;
	int k2=0;

	int start=-1;

	for(int i=0;i<length;i++){
		if(fabsf(freArr[i]-value)<10){
			start=i;
			break;
		}
	}

	if(start==-1||
		start+3>length-1){

		return 0;
	}

	flag=1;
	for(int i=start+1,j=2;i<length&&j<5;i++,j++){
		__queue_fre2(freArr[start], freArr[i], &k1, &k2);
		if(!(k1==1&&k2==j)){
			flag=0;
			break;
		}
	}

	return flag;
}

static int __queue_query2Inf(float *freArr,float *dbArr,float *heightArr,int length,int start,float value1,float value2,int strict){
	int flag=0;
	int count=0;

	for(int i=start;i<length;i++){
		int k1=0;
		int k2=0;
		int type=0;

		k1=util_calRangeTimes(value1, freArr[i], &type);
		if(k1&&!type){
			if(k1==3||k1==5||k1==7||k1==9||k1==11||k1==13){
				if(k1<=5&&fabsf(value1*k1-freArr[i])<6){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i-1]-dbArr[i]<3||
								dbArr[i+1]-dbArr[i]<3){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1<=9&&fabsf(value1*k1-freArr[i])<12){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i-1]-dbArr[i]<3||
								dbArr[i+1]-dbArr[i]<3){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1==11&&fabsf(value1*k1-freArr[i])<18){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i-1]-dbArr[i]<3||
								dbArr[i+1]-dbArr[i]<3){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1==13&&fabsf(value1*k1-freArr[i])<20){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i-1]-dbArr[i]<3||
								dbArr[i+1]-dbArr[i]<3){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
			}

			if(strict<2&&count){
				flag=1;
				break;
			}
			else if(strict>=2&&count>=strict){
				flag=1;
				break;
			}
		}
	}

	return flag;
}

static int __queue_valid98(float *freArr,float *dbArr,float *heightArr,int length,int start,float value1,int strict){
	int count=0;

	for(int i=start;i<length;i++){
		int k1=0;
		int k2=0;
		int type=0;

		k1=util_calRangeTimes(value1, freArr[i], &type);
		if(k1&&!type){
			if(k1==3&&i-1>=0){
				if(dbArr[i-1]-dbArr[i]>24){
					continue;
				}
			}

			if(k1==3||k1==5||k1==7||k1==9||k1==11||k1==13||k1==15||k1==17||k1==19){
				if(k1<=5&&fabsf(value1*k1-freArr[i])<6){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i]-dbArr[i-1]>2||
								dbArr[i]-dbArr[i+1]>2){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1<=9&&fabsf(value1*k1-freArr[i])<18){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i]-dbArr[i-1]>2||
								dbArr[i]-dbArr[i+1]>2){

								if(dbArr[i-1]-dbArr[i]>24){
									if(dbArr[i]-dbArr[i+1]>6){
										count++;
									}
								}
								else{
									count++;
								}
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1==11&&fabsf(value1*k1-freArr[i])<20){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i]-dbArr[i-1]>2||
								dbArr[i]-dbArr[i+1]>2){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1<=19&&fabsf(value1*k1-freArr[i])<25){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<3){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i]-dbArr[i-1]>3||
								dbArr[i]-dbArr[i+1]>3){

								if(dbArr[i-1]-dbArr[i]>18){
									if(dbArr[i]-dbArr[i+1]>6){
										count++;
									}
								}
								else{
									count++;
								}
							}
						}
					}
					else{
						count++;
					}
				}
			}
		}
	}

	return count;
}

static int __queue_odd98(float *freArr,float *dbArr,float *heightArr,int length,int start,float value1){
	int count=0;
	int corrFlag=0;

	for(int i=start;i<length;i++){
		int k1=0;
		int k2=0;
		int type=0;

		k1=util_calRangeTimes(value1, freArr[i], &type);
		if(k1&&!type){
			if(!corrFlag){
				if(k1==4||k1==6||k1==8){
					if(fabsf(value1*2-196)>fabsf(freArr[i]/k1*2-196)){
						value1=freArr[i]/k1;
					}

					corrFlag=1;
				}
			}
			
			if(k1%2==1&&k1>1){
				if(k1<=5&&fabsf(value1*k1-freArr[i])<6){
					count++;
				}
				else if(k1<=9&&fabsf(value1*k1-freArr[i])<18){
					count++;
				}
				else if(k1==11&&fabsf(value1*k1-freArr[i])<20){
					count++;
				}
				else if(k1<=19&&fabsf(value1*k1-freArr[i])<25){
					count++;
				}
				else if(k1>19&&fabsf(value1*k1-freArr[i])<30){
					count++;
				}
			}
		}
	}

	return count;
}

static int __queue_query2(float *freArr,float *dbArr,float *heightArr,int length,int start,float value1,float value2,int strict){
	int flag=0;
	int count=0;

	for(int i=start;i<length;i++){
		int k1=0;
		int k2=0;
		int type=0;

		k1=util_calRangeTimes(value1, freArr[i], &type);
		if(k1&&!type){
			if(k1==3||k1==5||k1==7||k1==9||k1==11||k1==13){
				if(k1<=5&&fabsf(value1*k1-freArr[i])<6){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<8){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<8){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i-1]-dbArr[i]<8||
								dbArr[i+1]-dbArr[i]<8){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1<=9&&fabsf(value1*k1-freArr[i])<12){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<8){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<8){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i-1]-dbArr[i]<8||
								dbArr[i+1]-dbArr[i]<8){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1==11&&fabsf(value1*k1-freArr[i])<18){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<8){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<8){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i-1]-dbArr[i]<8||
								dbArr[i+1]-dbArr[i]<8){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
				else if(k1==13&&fabsf(value1*k1-freArr[i])<20){
					if(strict){
						if(i==start&&i<length-1){
							if(dbArr[i+1]-dbArr[i]<8){
								count++;
							}
						}
						else if(i==length-1&&i>0){
							if(dbArr[i-1]-dbArr[i]<8){
								count++;
							}
						}
						else if(i>0&&i<length-1){
							if(dbArr[i-1]-dbArr[i]<8||
								dbArr[i+1]-dbArr[i]<8){

								count++;
							}
						}
					}
					else{
						count++;
					}
				}
			}

			if(strict<2&&count){
				flag=1;
				break;
			}
			else if(strict>=2&&count>=strict){
				flag=1;
				break;
			}
		}
	}

	return flag;
}

static int __queue_query3(float *freArr,float *dbArr,float *heightArr,int length,int start,float value1,float value2,int strict){
	int flag=0;

	for(int i=start;i<length;i++){
		int k1=0;
		int k2=0;
		int type=0;

		k1=util_calRangeTimes(value1, freArr[i], &type);
		if(k1&&!type){
			if(k1==4||k1==5||k1==7||k1==8||k1==10||k1==11||k1==13){
				if(k1<=5&&fabsf(value1*k1-freArr[i])<6){
					flag=1;
				}
				else if(k1<=9&&fabsf(value1*k1-freArr[i])<12){
					flag=1;
				}
				else if(k1<=11&&fabsf(value1*k1-freArr[i])<18){
					flag=1;
					k2=util_calRangeTimes(value2, freArr[i], &type);
					if(k1==10&&k2==3&&
						fabsf(value2*k2-freArr[i])<8){
						
						flag=0;
					}
				}
				else if(k1==13&&fabsf(value1*k1-freArr[i])<20){
					flag=1;
					k2=util_calRangeTimes(value2, freArr[i], &type);
					if(k2==4&&
						fabsf(value2*k2-freArr[i])<10){

						flag=0;
					}
				}
			}

			if(flag){
				break;
			}
		}
	}

	return flag;
}

/***
	1. 123-4->2468; 123-6->246-12; 1246->248-12
	2. 2x3/12x3 ->456/2456, x-min ???
	2.1 247->
	3. 1246/2346 ???

****/
static float __queue_cutValid(float *freArr,float *dbArr,int length,int oFlag,int mode,
							float *freArr2,float *dbArr2,int length2,
							float *freArr3,float *dbArr3,int length3){
	float fre=0;
	int count=0;

	float *_freArr=NULL;
	float *_dbArr=NULL;
	int _length=0;

	if(mode==0){
		_freArr=freArr2;
		_dbArr=dbArr2;
		_length=length2;
	}
	else{
		_freArr=freArr3;
		_dbArr=dbArr3;
		_length=length3;
	}

	fre=freArr[0];
	for(int i=0;i<_length;i++){
		int k=0;
		int type=0;

		k=util_calRangeTimes(freArr[0]/2, _freArr[i], &type);

		if(oFlag&&!mode){ // ->247
			if(k==1&&length3<6){
				count++;
			}
		}

		if(k&&!type){
			if((!mode&&(k==3||k==5||k==7))||
				(mode&&(k==3||k==5||k==7||k==9||k==11))){ // ||k==9||k==11

				int flag=0;

				if(k<=5&&fabsf(freArr[0]/2*k-_freArr[i])<6){
					flag=1;
				}
				else if(k<=9&&fabsf(freArr[0]/2*k-_freArr[i])<12){
					flag=1;
				}
				else if(k==11&&fabsf(freArr[0]/2*k-_freArr[i])<18){
					flag=1;
				}

				if(flag&&
					(_dbArr[i-1]-_dbArr[i]<12||
					_dbArr[i+1]-_dbArr[i]<12)){

					count++;
				}
			}
		}
	}

	if(count==1&&oFlag&&length3>5){
		int us1=0,us2=0;
		int uk1=0,uk2=0,uk3=0;

		float _fre=0;
		int _k=0;

		for(int i=3;i<_length-2;i++){
			if(i>5){
				break;
			}

			_fre=__queue_fre3(_freArr[i],_freArr[i+1],_freArr[i+2],
						&us1,&us2,
						&uk1,&uk2,&uk3);
			
			if(us1==1&&us1==us2){ // 1:1
				_k=util_calRangeTimes(_fre, freArr[0], NULL);
				if(_k==2&&fabsf(_fre-freArr[0]/2)<8){
					count++;
					break;
				}
			}
		}
	}

	if(count>=2){
		fre=freArr[0]/2;
	}
	
	return fre;
}

/***
	6string, 82.4 ->75~90, 150~180
	5string, 110 ->103~118
	4string, 147 ->140~154
	3string, 196 ->190~214
	2string, 247 ->240~254
	1string, 329.6 ->320~340
****/
float __queue_standard(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,
					float *freArr2,float *dbArr2,float *heightArr2,int length2,
					float *freArr3,float *dbArr3,float *heightArr3,int refLength,
					float light,int *valid,
					int *formatFlag,
					float *fre11,float *fre22,float *fre33,
					float *db11,float *db22,float *db33){

	float fre=0;
	int status=0;

	float *_corrArr2=NULL;
	float *_dbArr2=NULL;

	int *_indexArr2=NULL;

	if(refLength>3){
		_corrArr2=__vnew(refLength, NULL);
		_dbArr2=__vnew(refLength, NULL);

		__varangei(0, refLength, 1, &_indexArr2);

		memcpy(_corrArr2, freArr3,sizeof(float )*refLength);
		memcpy(_dbArr2, dbArr3,sizeof(float )*refLength);

		// dB desc ->corrArr2
		__vcorrsort1(_dbArr2, _corrArr2 ,NULL,NULL,refLength, 1);
		// fre asc ->indexArr1
		__vcorrsort1(_corrArr2, NULL ,NULL,_indexArr2,refLength, 0);

		fre=__queue_slide(freArr3,dbArr3,heightArr3,_indexArr2,refLength,light,valid,&status);

		if(fre>240){
			
		}
		else if(fre>230&&refLength>12){ // 2-string 230~240&&refLength>12
			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, fre/2, fre,0);
			if(flag){
				fre=0;
			}
		}
		else{
			if(fre>189&&fre<205&&
				refLength>13){ // 3-string 197+7&&refLength>13
				
				int flag=0;

				int k1=0;
				int index1=0;

				for(int i=0;i<refLength;i++){
					if(fabs(fre*2-freArr3[i])<10){
						k1=2;
						index1=i;

						break;
					}
					else if(fabs(fre*3-freArr3[i])<15){
						k1=3;
						index1=i;

						break;
					}
				}

				if(k1){
					flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr3[index1]/k1, freArr3[index1]/k1*2,0);
					if(flag){
						int c1=0;
						int c2=0;

						int count1=0;

						count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr3[index1]/k1/2);
						if(count1>3){
							fre=freArr3[index1]/k1/2;
						}
						else{
							// flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr3[index1]/k1/2, freArr3[index1]/k1,1);
							// flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr3[index1]/k1/2,1);

							// if(flag>1){
							// 	c1=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr3[index1]/k1/2*13+20, freArr3[index1]/k1/2, 0);
							// 	c2=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr3[index1]/k1/2*13+20, freArr3[index1]/k1, 1);
							// 	if(!c1&&c2>=1){
							// 		flag=0;
							// 	}

							// 	if(flag&&!c1&&
							// 		freArr[1]>190&&freArr[1]<205&&
							// 		dbArr[1]-dbArr[0]>24){

							// 		flag=0;
							// 	}

							// 	if(flag&&freArr[3]>2000){
							// 		flag=0;
							// 	}
							// }

							// if(flag>1){
							// 	fre=freArr3[index1]/k1/2;
							// }
							// else{
							// 	fre=freArr3[index1]/k1;
							// }

							fre=freArr3[index1]/k1;
						}
					}
					else{
						fre=0;
					}
				}
			}
			else if(fre>139&&fre<155&&
				refLength>15){ // 4-string 147

				int flag=0;

				flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, fre, fre*2,0);
				if(!flag){
					fre=0;
				}
			}
			else{
				fre=0;
			}
		}

		if(fre>280&&fre<310){ // valid ->147
			int flag=0;
			int count=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, fre/2, fre,0);
			if(flag){
				flag=1;
				if(freArr[0]>190&&freArr[0]<205){
					count=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, 0, fre, 2);
					if(count>=2){
						flag=0;
					}
				}

				if(flag){
					fre=fre/2;
				}
				else{
					fre=0;
				}
			}
			else{
				fre=0;
			}
		}

		if(((fre/2>190&&fre/2<205)||
				(fre/4>190&&fre/4<205))&&
			refLength>4){ // valid ->197, 391/782

			int index1=0;
			int index2=0;

			index1=__arr_maxIndex(dbArr3, refLength);
			index2=__arr_maxIndex(dbArr3+1, refLength-1)+1;
			if((freArr3[index1]>179&&freArr3[index1]<205)||
				(freArr3[index2]>179&&freArr3[index2]<205)){

				if(fre/2>190&&fre/2<205){
					fre=fre/2;
				}
				else{
					fre=fre/4;
				}
			}
		}

		if(((fre/2>240&&fre/2<255)||
				(fre/4>240&&fre/4<255))&&
			refLength>8){ // valid ->247, 492/984

			float _fre1=0;

			int _num=2;
			int _subType=0;
			int _unionType=2;
			int _direction=0;

			_fre1=__queue_multi(freArr3,dbArr3,heightArr3,refLength,_num,_subType, _unionType, _direction);
			if(_fre1>240&&_fre1<255){
				fre=_fre1;
			}
		}

		if(fre>310&&fre<350&&
			freArr[0]>100&&freArr[0]<120&&
			dbArr[0]-dbArr[2]>10){

			int flag=0;

			int us1=0,us2=0;
			int uk1=0,uk2=0,uk3=0;

			int k1=0;
			int k2=0;

			flag=__queue_query3(freArr3,dbArr3,heightArr3,refLength,0,freArr[0],fre,0);
			if(flag){ // ->110
				fre=freArr[0];
			}

			if(!flag){
				__queue_fre3(freArr[0],freArr[1],freArr[2],
							&us1,&us2,
							&uk1,&uk2,&uk3);

				if(uk1==1&&uk2==2&&uk3==3&&
					fabsf(freArr[0]*2-freArr[1])<5&&
					fabsf(freArr[0]*3-freArr[2])<5){

					fre=freArr[0];
				}
				else{
					__queue_fre2(freArr[1],freArr[2],
							&k1,&k2);

					if(k1==2&&k2==3&&
						fabsf(freArr[1]/2*3-freArr[2])<6&&
						fabsf(freArr[0]-freArr[1]/2)<8){

						fre=freArr[0];
					}
				}
			}
		}

		// if(fre>310&&fre<350&&
		// 	freArr[0]>100&&freArr[0]<120){

		// 	int flag=0;

		// 	int k1=0;
		// 	int k2=0;

		// 	__queue_fre2(freArr[0],fre,
		// 				&k1,&k2);

		// 	if(k1==1&&k2==3&&
		// 		fabsf(freArr[0]*3-fre)<4){

		// 		flag=__queue_query3(freArr3,dbArr3,heightArr3,refLength,0,freArr[0],fre);
		// 		if(flag){
		// 			fre=fre/3;
		// 		}
		// 	}
		// }

		free(_corrArr2);
		free(_dbArr2);
		free(_indexArr2);
	}
	
	return fre;
}

/***
	valid
		1:1 ->valid 1:2:4:6/2:3:4:6 & valid 1:x:2 & valid 2:x:3 
		456/2456 ->2x3,12x3
	normal
		1234/2345/3456
		2347/2467
		23+67/23+56/12+67/12+56

		367-11,
	area
		2346/1346(146/346)
		1345/1356/1357/1456/1457

	1234->2468; 1236->246-12; 1246->248-12
	236+7/8/10/11; max1-2

	1245/1345 ???
	349??? + valid ->1/5/7
	234-12 -> 1x26
****/
float __queue_cut(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,
				float *freArr2,float *dbArr2,float *heightArr2,int length2,
				float *freArr3,float *dbArr3,float *heightArr3,int refLength,
				float light,int *valid,
				int *formatFlag,
				float *fre11,float *fre22,float *fre33,
				float *db11,float *db22,float *db33){

	float fre=0;
	float _fre=0;
	
	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int vs1=0,vs2=0;
	int vk1=0,vk2=0,vk3=0;

	int index1=0;

	index1=__arr_maxIndex(dbArr, length);

	// 1x23 ->valid 234-6
	if((fabsf(dbArr[2]-dbArr[3])<4||
		dbArr[2]>dbArr[3])&&
			dbArr[2]>dbArr[0]&&dbArr[2]>dbArr[1]&&
			dbArr[3]>dbArr[0]&&dbArr[3]>dbArr[1]){ // ->110, ->82, ->70
		
		_fre=__queue_fre3(freArr[0],freArr[2],freArr[3],
						&us1,&us2,
						&uk1,&uk2,&uk3);

		__queue_fre3(freArr[0],freArr[1],freArr[2],
						&vs1,&vs2,
						&vk1,&vk2,&vk3);

		if(uk1==1&&uk2==2&&uk3==3){
			if(fabsf(_fre*uk2-freArr[2])<5&&
				fabsf(_fre*uk3-freArr[3])<5){

				if(vk2!=3){
					// fre=__checkFre(freArr[0]/uk1, freArr[2]/uk2, freArr[3]/uk3);
					fre=freArr[2]/uk2;
					return fre;
				}
				else{
					if(freArr[0]<100){
						// fre=__checkFre(freArr[0]/uk1, freArr[2]/uk2, freArr[3]/uk3);
						fre=freArr[2]/uk2;
						return fre;
					}
				}
			}
		}
	}
	else if(dbArr[0]-dbArr[1]>12&&
			dbArr[2]-dbArr[1]>12&&
			freArr[1]>160&&freArr[1]<180){ // 1x23 ->valid-110

		_fre=__queue_fre3(freArr[0],freArr[2],freArr[3],
						&us1,&us2,
						&uk1,&uk2,&uk3);

		if(uk1==1&&uk2==2&&uk3==3){
			if(fabsf(_fre*uk2-freArr[2])<5&&
				fabsf(_fre*uk3-freArr[3])<5){

				// fre=__checkFre(freArr[0]/uk1, freArr[2]/uk2, freArr[3]/uk3);
				fre=freArr[2]/uk2;
				return fre;
			}
		}
	}
	else if(freArr[0]>103&&freArr[0]<115){
		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if(!_index){
			_fre=__queue_fre3(freArr[0],freArr[2],freArr[3],
							&us1,&us2,
							&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==3&&uk3==4){ // 1x34 ->valid-110
				if(fabsf(_fre*uk2-freArr[2])<5&&
					fabsf(_fre*uk3-freArr[3])<5){

					if((freArr[0]*2-freArr[1])<15){
						// fre=__checkFre(freArr[0]/uk1, freArr[2]/uk2, freArr[3]/uk3);
						fre=freArr[2]/uk2;
						return fre;
					}
				}
			}
			else if(uk1==1&&uk2==4&&uk3==6){ // 1x46 ->valid-110
				if(fabsf(_fre*uk2-freArr[2])<5&&
					fabsf(_fre*uk3-freArr[3])<5){

					if((freArr[0]*2-freArr[1])<15){
						// fre=__checkFre(freArr[0]/uk1, freArr[2]/uk2, freArr[3]/uk3);
						fre=freArr[2]/uk2;
						return fre;
					}
				}
			}
		}
	}

	// 1234/2345/3456 & 23+67/23+56/12+67/12+56
	__queue_fre3(freArr[0],freArr[1],freArr[2],
				&us1,&us2,
				&uk1,&uk2,&uk3);
	if(uk1){
		__queue_fre3(freArr[1],freArr[2],freArr[3],
					&vs1,&vs2,
					&vk1,&vk2,&vk3);
		if(vk1){
			if(uk1>=1&uk1<=2&uk1+1==uk2&&uk2+1==uk3&&uk3+1==vk3){ // ->2345
				// fre=freArr[0]/uk1;
				// fre=__checkFre(freArr[0]/uk1, freArr[1]/uk2, freArr[2]/uk3);
				fre=freArr[1]/uk2;

				if(uk1==1){ // ->1234 cut_valid
					if(dbArr[0]>dbArr[1]&&
						(dbArr[1]>dbArr[2]&&
							dbArr[1]>dbArr[3])){

						float _fre1=0;
						int k1=0,k2=0;

						_fre1=__queue_cutValid(freArr,dbArr,length,0,1,
											freArr2,dbArr2,length2,
											freArr3,dbArr3,refLength);

						__queue_fre2(_fre1, fre, &k1, &k2);
						if(!(k1==1&&k1==k2)){
							fre=_fre1;
						}
					}
					else{ 
						// if(dbArr[2]>dbArr[1]&&
						// 	dbArr[1]>dbArr[0]&&
						// 	dbArr[0]>dbArr[3]&&
						// 	freArr[0]>120&&freArr[0]<160){ // area ->74

						// 	fre=freArr[0]/2;
						// }

						if(index1==1&&
							freArr[index1]>190&&freArr[index1]<204&&
							dbArr[2]<dbArr[0]&&
							heightArr[2]<15){ // 234 ->197

							return freArr[1];
						}
						else if(index1==1&&
							freArr[index1]>190&&freArr[index1]<204&&
							dbArr[1]-dbArr[2]>18){ // 234 ->197

							int count1=0;

							count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
							if((count1>=2&&refLength<10)||
								count1>=3){

								return freArr[1]/2;
							}

							if(refLength<7){
								return freArr[1]/2;
							}
							else{
								int flag=0;

								int count1=0;
								int flag1=0;

								flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, 1);
								count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
								if(count1>=2&&refLength<10){
									flag1=1;
								}
								else if(count1>3){
									flag1=1;
								}
								else if(dbArr[1]-dbArr[0]<6&&
									refLength<10&&
									count1){

									flag1=1;
								}

								if(!flag&&!flag1){
									return freArr[1];
								}
							}
						}
						else if(freArr[0]>150&&freArr[0]<180){ // 2468 ->80{75,90}
							int flag=0;

							// query2
							flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[0]/2, freArr[0],1);
							if(flag){
								return freArr[0]/2;
							}
						}
					}
				}

				return fre;
			}

			if(uk1==2&&uk2==3&&uk3==4&&vk3==7){ // ->2347
				fre=freArr[0]/uk1;
				return fre;
			}

			if(uk1==1&&uk2==2&&uk3==3&&
				vk1==4&&vk2==6&&vk3==7){ // ->2467

				if(dbArr[0]>dbArr[1]&&
					dbArr[1]>dbArr[2]&&
					dbArr[2]>dbArr[3]&&
					freArr[0]>100&&freArr[0]<120){ // valid ->110

					fre=freArr[1]/2;
				}
				else if(index1==1&&
					dbArr[index1]-dbArr[3]>20&&
					freArr[1]<190){
					
					fre=freArr[1]/2;
				}
				else if((index1==1||index1==2)&&
					dbArr[1]-dbArr[3]>12&&
					dbArr[2]-dbArr[3]>12&&
					freArr[1]>150&&freArr[1]<180){ // valid ->80, 2467,123x

					fre=freArr[1]/2;
				}
				else{
					if(dbArr[0]-dbArr[3]>20&&
						(dbArr[1]-dbArr[3]>20||
							dbArr[0]-dbArr[1]>20)){ // valid ->246(7), 123x

						fre=freArr[0];

					}
					else if(dbArr[0]-dbArr[3]>18&&
						dbArr[1]-dbArr[3]>12&&
						dbArr[2]>dbArr[3]&&
						freArr[0]>220){ // ->247,123x

						fre=freArr[1]/2;
					}
					else{
						fre=freArr[0]/2;
					}
				}
				
				return fre;
			}

			if((uk1==2&&uk2==3&&uk3==6&&
					vk1==3&&vk2==6&&vk3==7)||
				(uk1==2&&uk2==3&&uk3==5&&
					vk1==3&&vk2==5&&vk3==6)||
				(uk1==1&&uk2==2&&uk3==6&&
					vk1==2&&vk2==6&&vk3==7)||
				(uk1==1&&uk2==2&&uk3==5&&
					vk1==2&&vk2==5&&vk3==6)){ // ->23+67/23+56/12+67/12+56

				if((uk1==2&&uk2==3&&uk3==5&&
					vk1==3&&vk2==5&&vk3==6)&&
					dbArr[1]>dbArr[3]&&
					dbArr[3]>dbArr[0]&&
					dbArr[3]>dbArr[2]&&
					freArr[1]>140&&freArr[1]<155){ // ->147, ❌48

					fre=freArr[1];
				}
				else{
					fre=freArr[0]/uk1;
				}
				
				return fre;
			}
		}
	}

	if(uk1==1&&uk2==2&&uk3==4&&
		vk1==1&&vk2==2&&vk3==4&&
		index1==2&&
		freArr[1]>103&&freArr[1]<120){ // ->110, n124,max-2

		return freArr[2]/2;
	}

	// x234 ->110, YAMAFA-110!!!, 23xx ->80, 12/1n ->329, 1x23 ->147, x234 ->80, x123 ->246, x23n ->196
	if(!uk1){
		int ts1=0,ts2=0;
		int tk1=0,tk2=0,tk3=0;

		int k1=0,k2=0;

		__queue_fre3(freArr[1],freArr[2],freArr[3],
						&ts1,&ts2,
						&tk1,&tk2,&tk3);

		if(tk1==2&&tk2==3&&tk3==4&&
			freArr[1]/2>100&&freArr[1]/2<120){ // ->110

			return freArr[1]/2;
		}

		if(tk1==1&&tk2==2&&tk3==3&&
			freArr[1]/2>100&&freArr[1]/2<120){ // ->110

			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1], 0);
			if(flag){
				return freArr[1]/2;
			}
		}

		if(tk1==2&&tk2==4&&tk3==7&&
			freArr[1]/2>100&&freArr[1]/2<120){ // ->110

			return freArr[1]/2;
		}

		if(tk1==2&&tk2==3&&tk3==4&&
			index1==1&&
			freArr[0]>85&&freArr[0]<95&&
			freArr[1]>150&&freArr[1]<170){ // x234 ->80, 75~85

			return freArr[1]/2;
		}

		if(index1==2&&
			dbArr[2]-dbArr[1]>18){ // ->147

			__queue_fre3(freArr[0],freArr[2],freArr[3],
						&ts1,&ts2,
						&tk1,&tk2,&tk3);

			if(tk1==1&&tk2==2&&tk3==3&&
				freArr[0]>140&&freArr[0]<154){

				return freArr[2]/2;
			}

			if(tk1==1&&tk2==3&&tk3==4&&
				freArr[2]>200&&freArr[0]<210){ // ->71 ,67~71

				return freArr[0];
			}
		}

		if(tk1==1&&tk2==2&&tk3==4&&
			index1==3&&
			freArr[1]>220&&freArr[1]<360){ // ->330, 124

			return freArr[2]/2;
		}

		if(tk1==1&&tk2==2&&tk3==4&&
			index1==2&&
			freArr[2]/2>140&&freArr[2]/2<155){ // ->147, 124

			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2, freArr[2], 0);
			if(flag){
				return freArr[2]/2;
			}
		}

		if(tk1==1&&tk2==2&&tk3==4&&
			index1==2&&
			freArr[2]/2>105&&freArr[2]/2<115){ // ->110, 124

			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2, freArr[2], 0);
			if(flag){
				return freArr[2]/2;
			}
		}

		__queue_fre2(freArr[1],freArr[2],
					&k1,&k2);

		if(index1==2&&
			dbArr[2]-dbArr[1]>18&&
			freArr[2]>300&&freArr[2]<350){ // ->330

			int flag=0;

			int _k1=0;
			int _k2=0;

			__queue_fre2(freArr[0],freArr[2],
						&_k1,&_k2);

			if(_k1==1&&_k2==3&&
				fabsf(freArr[0]*3-freArr[2])<4){

				flag=__queue_query3(freArr3,dbArr3,heightArr3,refLength,0,freArr[0],freArr[2],0);
				if(flag){
					return freArr[2]/3;
				}
			}

			return freArr[2];
		}

		if(k1==2&&k2==3&&
			freArr[1]>200&&freArr[1]<240&&
			fabsf(freArr[1]/2*3-freArr[2])<4&&
			dbArr[1]-dbArr[2]>-10){

			if(fabsf(freArr[0]-freArr[1]/2)<15||
				(freArr[0]>90&&freArr[0]<100)){ // (freArr[0]>90&&freArr[0]<110)

				return freArr[1]/2;
			}
		}

		if(index1==0&&
			dbArr[1]>dbArr[3]&&
			dbArr[2]>dbArr[3]&&
			freArr[2]/3>100&&freArr[2]/3<120){ // 123 ->1x3

			__queue_fre2(freArr[0],freArr[1],
						&k1,&k2);

			if(k1==1&&k2==2&&
				fabsf(freArr[1]/2-freArr[0])<6){

				int tk1=0,tk2=0;

				__queue_fre2(freArr[0],freArr[2],
							&tk1,&tk2);

				if(tk1==1&&tk2==3){
					return freArr[2]/3;
				}
			}
		}

		if(index1==0&&
			dbArr[2]>dbArr[1]&&
			dbArr[3]>dbArr[1]){ // ->110

			__queue_fre3(freArr[0],freArr[2],freArr[3],
						&ts1,&ts2,
						&tk1,&tk2,&tk3);

			if(tk1==1&&tk2==2&&tk3==3&&
				freArr[2]>200&&freArr[2]<240){

				return freArr[0];
			}
		}

		if((index1==1||index1==0)&&
			fabsf(dbArr[0]-dbArr[1])<3&&
			dbArr[0]>dbArr[2]&&
			dbArr[1]>dbArr[2]){ // ->110, deform{118,229,330,456}

			if(freArr[0]>110&&freArr[0]<120&&
				freArr[1]>220&&freArr[1]<240&&
				freArr[2]>315&&freArr[2]<345&&
				freArr[3]>420&&freArr[3]<460){

				return freArr[2]/3;
			}
		}

		if(index1==1&&
			tk1==1&&tk2==2&&tk3==3&&
			freArr[2]/2>230&&freArr[2]/2<255){ // ->246, x123

			return freArr[2]/2;
		}

		if(index1==2&&
			tk1==1&&tk2==2&&(tk3==4||tk3==6)&&
			freArr[2]/2>95&&freArr[2]/2<105){ // ->100, x124/x126

			int flag=0;
			int count1=0;

			// flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2, freArr[2],1);
			flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2, 1);
			if(!flag&&refLength<8){
				count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2);
				if(count1>=2){
					flag=1;
				}
			}

			if(flag){
				return freArr[2]/2;
			}
			else{
				return freArr[2];
			}
		}

		if(index1==1&&
			tk1==2&&tk2==3&&tk3==6&&
			freArr[1]/2>95&&freArr[1]/2<105){ // ->100, x236 ❌❌❌

			int flag=0;
			int count1=0;

			// flag=__queue_query2(freArr3+_len, dbArr3+_len, heightArr3+_len, refLength-_len, 0, freArr[1]/2, freArr[1],1);
			flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2,1);
			count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
			if(!flag&&refLength<10){
				count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
				if(count1>=2){
					flag=1;
				}
			}

			if(flag||count1>3){
				return freArr[1]/2;
			}
			else{
				return freArr[1];
			}
		}

		if(index1==1&&
			tk1==2&&tk2==3&&tk3==4&&
			freArr[1]/2>95&&freArr[1]/2<105&&
			freArr[0]<100){ // ->100, x234

			int flag=0;
			int count1=0;

			flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2,1);
			if(!flag&&refLength<8){
				count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
				if(count1>=2){
					flag=1;
				}
			}

			if(flag){
				return freArr[1]/2;
			}
			else{
				return freArr[3]/2;
			}
		}

		if(index1==1&&
			tk1==1&&tk2==2&&tk3==3&&
			freArr[1]/2>95&&freArr[1]/2<105&&
			freArr[0]<110){ // ->100, x246

			int flag=0;
			int count1=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],1);
			if(flag){
				int c1=0;
				int c2=0;

				c1=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2*11+10, freArr[1]/2, 0);
				c2=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2*11+10, freArr[1], 2);
				if(!c1&&c2>=1){
					flag=0;
				}

				if(flag){
					flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2,1);
				}
			}

			if(!flag&&refLength<8){
				count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
				if(count1>=2){
					flag=1;
				}
			}

			if(flag){
				return freArr[1]/2;
			}
			else{
				return freArr[2]/2;
			}
		}

		if(index1==2&&
			!tk1&&
			freArr[2]/2>95&&freArr[2]/2<105&&
			freArr[1]>95&&freArr[1]<106){ // ->100, xx12

			int flag=0;
			int count1=0;

			int _k1=0;
			int _k2=0;

			__queue_fre2(freArr[2],freArr[3],
						&_k1,&_k2);

			if(_k1==1&&_k2==2&&
				fabsf(freArr[2]*2-freArr[3])<5){

				flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2, freArr[2],1);
				if(!flag&&refLength<8){
					count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2);
					if(count1>=2){
						flag=1;
					}
				}

				if(flag){
					return freArr[2]/2;
				}
				else{
					return freArr[3]/2;
				}
			}
		}

		if(index1==1&&
			!tk1&&
			freArr[2]>freArr[0]*6&&
			freArr[1]/2>95&&freArr[1]/2<105&&
			freArr[0]>92&&freArr[0]<106){ // ->100, 12nn

			int flag=0;
			int count1=0;

			int _k1=0;
			int _k2=0;

			__queue_fre2(freArr[0],freArr[1],
						&_k1,&_k2);

			if(_k1==1&&_k2==2){

				// flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],1);
				flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, 1);
				if(!flag&&refLength<8){
					__queue_fre2(freArr[0], freArr[2], &_k1, &_k2);
					if(_k1==1){
						float _fre=0;

						if(fabsf(freArr[1]-196)<fabsf(freArr[2]/_k2*2-196)){
							_fre=freArr[1]/2;
						}
						else{
							_fre=freArr[2]/_k2;
						}

						count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, _fre);
						if(count1>=2){
							flag=1;
						}
					}
				}

				if(flag){
					return freArr[1]/2;
				}
				else{
					return freArr[1];
				}
			}
		}

		if(index1==1&&
			!tk1&&
			freArr[1]/2>95&&freArr[1]/2<105&&
			freArr[0]>95&&freArr[0]<106){ // ->100, x13n

			int flag=0;
			int count1=0;

			int _k1=0;
			int _k2=0;

			__queue_fre2(freArr[1],freArr[2],
						&_k1,&_k2);

			if(_k1==1&&_k2==3&&
				fabsf(freArr[1]*3-freArr[2])<8){

				// flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],1);
				flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, 1);
				if(!flag&&refLength<8){
					count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
					if(count1>=2){
						flag=1;
					}
				}

				if(flag){
					return freArr[1]/2;
				}
				else{
					return freArr[2]/3;
				}
			}
		}

		// if(index1==1&&
		// 	dbArr[0]>dbArr[3]&&
		// 	dbArr[2]>dbArr3[3]&&
		// 	freArr[1]>200&&freArr[1]<240){ // 123 ->12x

		// 	if(fabsf(freArr[1]/2-freArr[0])<5&&
		// 		fabsf(freArr[2]/3-freArr[1]/2)<5){

		// 		return freArr[2]/3;
		// 	}
		// }

		// if(k1==2&&k2==3&&
		// 	dbArr[0]>dbArr[1]&&
		// 	dbArr[1]>dbArr[2]&&
		// 	dbArr[1]>dbArr[3]&&
		// 	freArr[0]>140&&freArr[0]<180){ // ->80

		// 	return freArr[0]/2;
		// }
	}

	// !uk1&&!vk1
	if(!uk1&&!vk1&&
		(index1==0||index1==1)&&
		freArr[1]>179&&freArr[1]<190&&
		freArr[2]/2>190&&freArr[2]/2<205){ // ->197, taylor

		int k1=0;
		int k2=0;

		__queue_fre2(freArr[2], freArr[3], &k1, &k2);
		if(k1==1&&k2==2&&
			fabsf(freArr[2]*2-freArr[3])<5){

			return freArr[2]/2;
		}
		else if(k1==2&&k2==3&&
			fabsf(freArr[2]/2*3-freArr[3])<5){

			return freArr[2]/2;
		}
	}

 	// 1245 ->valid 197, ->110
	if(uk1==1&&uk2==2&&uk3==4&&
		vk1==2&&vk2==4&&vk3==5){

		if(index1==1&&
			freArr[index1]>190&&freArr[index1]<204){

			return freArr[0];
		}

		if(index1==2&&
			freArr[1]/2>105&&freArr[1]/2<115){

			return freArr[1]/2;
		}
	}

	// 124, ->196, 124n, 4-max, ->110, 1246
	if(uk1==1&&uk2==2&&uk3==4){
		if(index1==2&&
			freArr[0]>185&&freArr[0]<205){

			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],0);
			if(flag){
				return freArr[1]/2;
			}
		}

		if(vk3==3&&
			index1==1&&
			freArr[0]>94&&freArr[0]<120){ // ->110 ,1246; ->100 ,1246

			int flag=0;
			int _index=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],0);
			if(flag){
				if(freArr[1]>207&&freArr[1]<230){ // ->110, 1246
					return freArr[1]/2;
				}
				else{ // ->100
					int c1=0;
					int c2=0;

					int count1=0;

					count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
					if(count1>=3){
						return freArr[1]/2;
					}

					c1=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2*13+20, freArr[1]/2, 0);
					c2=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2*13+20, freArr[1], 1);
					if(!c1&&c2>=1){
						return freArr[2]/2;
					}

					flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2,1);
					if(!flag){
						return freArr[2]/2;
					}
				}
			}

			flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2,1);
			if(!flag){
				return freArr[2]/2;
			}

			flag=__queue_bear(freArr3, dbArr3, heightArr3, refLength,  freArr[1]/2*13, freArr[1]/2, &_index);
			if(flag){
				return freArr[1]/2;
			}
		}
	}

	// 124 ->valid 110
	if(uk1==1&&uk2==2&&uk3==4&&
		!(vk1==2&&vk2==4&&vk3==5)){

		if(index1==1&&
			dbArr[1]-dbArr[0]>15&&
			freArr[index1]>100&&freArr[index1]<120){

			return freArr[2]/2;
		}
	}

	// 123n ->100
	if(uk1==1&&uk2==2&&uk3==3&&
		((dbArr[1]-dbArr[0]<6&&
			dbArr[1]-dbArr[2]>8)||
		(dbArr[0]-dbArr[1]>5&&
			dbArr[1]>dbArr[2]))&&
		freArr[0]>95&&freArr[0]<105){ // ->100

		return freArr[1]/2;
	}

	// // 124n ->100
	// if(uk1==1&&uk2==2&&uk3==4&&
	// 	dbArr[1]-dbArr[0]<3&&
	// 	dbArr[1]-dbArr[2]>10&&
	// 	freArr[0]>95&&freArr[0]<105){ // ->100

	// 	return freArr[1]/2;
	// }

	// 124/126/136/123 ->valid 1:2,1:3, ->220, ->98, ->294, ->147 
	{
		// 124n/126n ->valid 1:2, 110, 98
		if(uk1==1&&uk2==2&&
			(uk3==4||uk3==6)&&
			index1==1){ // freArr[1]>210&&freArr[1]<230

			int flag=0;
			int cutLen=0;

			int count1=0;
			int flag1=0;

			float _fre=0;

			if(vk1==2&&vk2==4&&vk3==5&&
				fabsf(freArr[0]*2-freArr[1])<5&&
				freArr[0]<95&&
				dbArr[1]-dbArr[0]<12&&
				dbArr[0]>dbArr[2]&&
				dbArr[0]>dbArr[3]){ // 1245 ->70~95

				return freArr[1]/2;
			}

			if(vk1==1&&vk2==2&&
				dbArr[1]-dbArr[0]>24&&
				freArr[1]>190&&freArr[1]<205){ // 1246/124-12, max-1, ->196

				return freArr[2]/2;
			}

			// ->147, ->196
			if(freArr[1]/2>140&&freArr[1]/2<155){
				return freArr[1]/2;
			}
			else if(freArr[1]/2>190&&freArr[1]/2<205){
				return freArr[1]/2;
			}

			// ->110
			if(freArr[1]/2>105&&freArr[1]/2<115){
				flag=1;
			}
			else if(freArr[1]/2>240&&freArr[1]/2<255){
				flag=1;
			}

			cutLen=__arr_cut(freArr3, refLength, freArr[1]*6);

			_fre=freArr[1]/2;
			if(fabsf(freArr[0]*uk3-freArr[2])<fabsf(freArr[1]*uk3/2-freArr[2])){
				_fre=freArr[0];
			}

			if(freArr[1]>190&&freArr[1]<205){
				flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, _fre,1);
				count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, _fre);

				if(dbArr[1]-dbArr[2]>20){
					flag=0;
				}

				if(count1>=2&&
					refLength<8){
					
					flag1=1;
				}
				else if(count1>3){
					flag1=1;
				}
				else if(dbArr[1]-dbArr[0]<6&&
					refLength<10&&
					count1){

					flag1=1;
				}
			}
			else{
				flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, _fre, _fre*2,flag?0:1);
			}

			if(!flag&&
				cutLen<6&&
				freArr[1]/2>105&&freArr[1]/2<115){ // ->110,

				flag=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]*7, _fre, 0);
			}

			if(flag||flag1){
				return freArr[1]/2;
			}
			else{
				return freArr[2]/(uk3/2);
			}
		}

		// 136n ->valid 1:3, !98{95,105}, 147->294->98
		if(uk1==1&&uk2==3&&uk3==6&&
			freArr[0]>95&&freArr[0]<105){

			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],1);
			if(flag){ // ->147
				return freArr[1]/2;
			}
			else{
				flag=__queue_query3(freArr3, dbArr3, heightArr3, refLength, 0, freArr[0], freArr[1],1);
				if(flag){ // ->98
					return freArr[1]/3;
				}
				else{ // ->294
					return freArr[2]/2;
				}
			}
		}

		// 123n ->valid 1:3, !98{95,105},{93,103}, ->98,->147,->196,->294
		if(uk1==1&&uk2==2&&uk3==3&&
			freArr[2]>280&&freArr[2]<310){

			int flag=0;
			int count1=0;
			int flag1=0;

			count1=__queue_odd98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2);
			
			flag=__queue_query3(freArr3, dbArr3, heightArr3, refLength, 0, freArr[0], freArr[2],1);
			if(flag||index1==1){ // ->98/196
				// flag=__queue_query2(freArr3+_len, dbArr3+_len, heightArr3+_len, refLength-_len, 0, freArr[1]/2, freArr[1],1);
				flag=__queue_valid98(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, 1);
				
				if(count1>=2&&
					refLength<8){

					flag1=1;
				}
				else if(count1>=3){
					flag1=1;
				}
				else if(dbArr[1]-dbArr[0]<6&&
					refLength<10&&
					count1){

					flag1=1;
				}

				if(flag||
					flag1||
					(vk1==2&&vk2==3&&
						(vk3==5||vk3==7))||
					(index1==1&&
						dbArr[1]-dbArr[2]>18&&
							(dbArr[2]-dbArr[3]>2||
								(dbArr[2]>dbArr[3]&&
									fabsf(freArr3[2]-freArr[2])<10)))){ // ->98

					return freArr[1]/2;
				}
				else{ // ->196
					return freArr[1];
				}
			}
			else{ // ->294
				return freArr[2];
			}
		}

		// 236n ->valid 1:3
		if(uk1==2&&uk2==3&&uk3==6&&
			index1>=1&&
			dbArr[index1]-dbArr[1]<3&&
			dbArr[1]>dbArr[0]&&
			freArr[0]/2>95&&freArr[0]/2<105){

			int flag=0;

			// flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],1);
			// if(flag){ // ->147
			// 	return freArr[1]/2;
			// }
			// else{
			// 	flag=__queue_query3(freArr3, dbArr3, heightArr3, refLength, 0, freArr[0]/2, freArr[1],1);
			// 	if(flag){ // ->98
			// 		return freArr[0]/2;
			// 	}
			// 	else{ // ->294
			// 		return freArr[2]/2;
			// 	}
			// }

			flag=__queue_query3(freArr3, dbArr3, heightArr3, refLength, 0, freArr[0]/2, freArr[1],1);
			if(flag){ // ->98
				return freArr[0]/2;
			}
			else{ // ->294
				return freArr[2]/2;
			}
		}

		// !uk1&&!vk1 ->294, x23x
		if(!uk1&&!vk1&&
			index1==1&&
			freArr[1]/2>280&&freArr[1]/2<310){

			int k1=0;
			int k2=0;

			__queue_fre2(freArr[1], freArr[2], &k1, &k2);
			if(k1==2&&k2==3&&
				fabsf(freArr[1]/2*3-freArr[2])<5){

				return freArr[1]/2;
			}
		}

		// ->294, x236, 
		if(index1<=2&&
			freArr[index1]>280&&freArr[index1]<310&&
			freArr[index1+1]/2>280&&freArr[index1+1]/2<310){ // 

			if(index1==2&&
				freArr[1]>140&&freArr[1]<155){

			}
			else{
				int count=0;

				count=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, 0, freArr[index1], 2);
				if(count>=2){
					return freArr[index1+1]/2;
				}
			}
		}
	}

	// 1367 ->valid 110, ->82
	if(uk1==1&&uk2==3&&uk3==6&&
		vk1==3&&vk2==6&&vk3==7){

		if(freArr[0]>100&&freArr[0]<120){

		}

		if(freArr[0]>75&&freArr[0]<90){
			return freArr[1]/3;
		}
	}

	// 245x ->valid 197, 147,2457
	if(uk1==2&&uk2==4&&uk3==5){
		if(index1==0&&
			freArr[index1]>190&&freArr[index1]<204){

			int _flag=0;

			for(int i=0;i<refLength;i++){
				if(fabsf(freArr[2]-freArr3[i])<1){
					_flag=1;
					break;
				}
			}

			if(_flag){
				return freArr[0]/2;
			}
			else{
				return freArr[0];
			}

			// if(dbArr[0]-dbArr[2]>20&&heightArr[2]<20){
			// 	return freArr[0];
			// }
			// else{
			// 	return freArr[0]/2;
			// }
		}
		else if(freArr[0]>280&&
			freArr[0]<310){ // 147, 2457

			return freArr[0]/2;
		}
	}

	// 145-8, valid-> 80 {75,90},1457
	if(uk1==1&&uk2==4&&uk3==5){
		if(index1==0&&
			freArr[0]>150&&freArr[0]<180){

			int flag=0;

			// query2
			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[index1]/2, freArr[index1],1);
			if(flag){
				return freArr[index1]/2;
			}
		}

		if(dbArr[0]-dbArr[1]>15){
			return freArr[0];
		}
	}

	// 1246 ->valid-197
	if(uk1==1&&uk2==2&&uk3==4&&
		vk3==3){ 

		if(freArr[1]>190&&freArr[1]<204&&
			dbArr[0]-dbArr[1]<3){

			if(fabsf(freArr[0]*2-freArr[1])<5){
				// fre=__queue_cutValid(freArr+1,dbArr+1,length-1,0,1,
				// 					freArr2,dbArr2,length2,
				// 					freArr3,dbArr3,refLength);

				int flag=0;

				flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],1);
				if(flag){
					fre=freArr[1]/2;
				}
				else{
					fre=freArr[2]/2;
				}

				return fre;
			}
			else{
				if(dbArr[1]>dbArr[2]&&
					dbArr[1]>dbArr[3]){

					return freArr[1];
				}
				else if(dbArr[2]>dbArr[1]&&
					dbArr[2]>dbArr[3]){

					return freArr[1];
				}
			}
		}

		if(freArr[1]>190&&freArr[1]<204){
			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],2);
			if(flag){
				return freArr[1]/2;
			}
			else{
				return freArr[2]/2;
			}
		}
	}

	// 124-12 ->valid-197
	if(uk1==1&&uk2==2&&uk3==4&&
		vk3==6){ 

		if(freArr[1]>190&&freArr[1]<204){
			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],1);
			if(flag){
				return freArr[1]/2;
			}
			else{
				return freArr[2]/2;
			}
		}
	}

	// 146 ->196, x23
	if(uk1==1&&uk2==4&&uk3==6){ 

		if(freArr[1]/2>190&&freArr[1]/2<204){
			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[0], freArr[1]/2,1);
			if(flag){
				return freArr[0];
			}
			else{
				return freArr[1]/2;
			}
		}
	}

	// 1236/1246 ->cut_valid, valid 80{75,90}
	if((uk1==1&&uk2==2&&uk3==3&&
		vk3==6)||
		(uk1==1&&uk2==2&&uk3==4&&
			vk3==3)){

		if(freArr[0]>75&&freArr[0]<90&&
			uk3==3){ // ->82, 1236,3-max/6-max

			return freArr[1]/2;
		}

		if(freArr[0]>95&&freArr[0]<105&&
			uk3==3&&
			(index1==2||index1==3)){ // ->!98, 1236,xx12

			int flag=0;

			if(index1==2&&
				freArr[2]>280&&freArr[2]<310){ // ->147

				int flag=0;

				flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2, freArr[2],0);
				if(flag){
					return freArr[2]/2;
				}
			}

			flag=__queue_query3(freArr3, dbArr3, heightArr3, refLength, 0, freArr[0], freArr[2],1);
			if(flag){
				return freArr[1]/2;
			}
			else{
				return freArr[3]/2;
			}
		}

		if(index1==0&&
			freArr[0]>150&&freArr[0]<180){

			int flag=0;

			// query2
			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[index1]/2, freArr[index1],1);
			if(flag){
				return freArr[index1]/2;
			}
		}

		if(dbArr[0]>dbArr[1]&&
			(dbArr[1]>dbArr[2]&&
				dbArr[1]>dbArr[3])){

			fre=__queue_cutValid(freArr,dbArr,length,0,0,
								freArr2,dbArr2,length2,
								freArr3,dbArr3,refLength);

			return fre;
		}
		else{ // valid ->110, ->80, ->197, x2x6, ->330, 1236, ->147,1236,xx36
			int _index=0;

			_index=__arr_maxIndex(dbArr, length);

			if(_index==2&&uk3==3&&
				dbArr[3]>dbArr[0]&&
				dbArr[3]>dbArr[1]&&
				freArr[2]>190&&freArr[2]<204){ // ->197

				return freArr[2];
			}

			if((_index==1||_index==2)&&uk3==3&&
				freArr[1]>130&&freArr[1]<180){ // ->80

				return freArr[1]/2;
			}

			// if(freArr[1]>200&&freArr[1]<240&&
			// 	dbArr[0]>dbArr[1]&&
			// 	dbArr[2]>dbArr[1]&&
			// 	uk3==3){ // ->110

			// 	int flag=0;

			// 	flag=__queue_query3(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/3, freArr[2],0);
			// 	if(flag){
			// 		return freArr[2]/3;
			// 	}
			// 	else{
			// 		return freArr[2];
			// 	}
			// }

			if(_index==2&&
				dbArr[2]-dbArr[1]>18&&
				uk3==3){ // ->330

				int flag=0;

				flag=__queue_query3(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/3, freArr[2],0);
				if(flag){
					return freArr[2]/3;
				}
				else{
					return freArr[2];
				}
			}

			if(uk3==3&&
				index1==2&&
				freArr[2]>280&&freArr[2]<310){ // ->147

				int flag=0;

				flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2, freArr[2],0);
				if(flag){
					return freArr[2]/2;
				}
			}

			if(index1==1&&uk3==4){ // ->110, 1246,x246
				int flag=0;

				flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[1]/2, freArr[1],0);

				if(!flag&&
					dbArr[1]-dbArr[0]<2&&
					fabsf(freArr[0]*2-freArr[1])<2){

					return freArr[1]/2;
				}

				if(flag){
					return freArr[1]/2;
				}
				else{
					return freArr[2]/2;
				}
			}

			if((dbArr[0]>dbArr[1]||
				dbArr[1]-dbArr[0]<3)&&
				(uk3==3?dbArr[2]-dbArr[1]>-10:1)&&
				freArr[0]>100&&freArr[0]<120){ // ->110

				return freArr[1]/2;
			}

			if(freArr[0]>186&&freArr[0]<206&&
				dbArr[1]>dbArr[0]&&
				dbArr[1]>dbArr[2]&&
				uk3==3){ // ->196, 1236

				return freArr[1]/2;
			}
		}
	}

	// 247/2478 ->cut_valid
	if(uk1==2&&uk2==4&&uk3==7&&
		fabsf(freArr[0]/2*7-freArr[2])<10){

		if(dbArr[0]>dbArr[1]&&
			(dbArr[1]>dbArr[2]&&
				dbArr[1]>dbArr[3])&&
			dbArr[0]-dbArr[2]>20){ // valid 24(7),12x

			return freArr[0];
		}

		if(dbArr[0]>dbArr[1]&&
			(dbArr[1]>dbArr[2]&&
				dbArr[1]>dbArr[3])){

			fre=__queue_cutValid(freArr,dbArr,length,1,0,
								freArr2,dbArr2,length2,
								freArr3,dbArr3,refLength);

			return fre;
		}
	}

	// 236 ->max0/2times valid ->147, 2369, ->197, ->80, ->71, 1x3, ->247,1x3 ->329,x12, ->247,x12
	if(uk1==2&&uk2==3&&uk3==6){
		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if(freArr[1]>130&&freArr[1]<160){ // ->147
			return freArr[1];
		}

		if(_index==0&&
			(dbArr[2]>dbArr[1]||(dbArr[0]-dbArr[1]>14&&heightArr[1]<20))&&
			freArr[0]>190&&freArr[0]<204){ // ->197

			return freArr[0];
		}

		if(_index==1&&
			freArr[1]>190&&freArr[1]<204){ // ->197, 130~140format

			return freArr[2]/2;
		}

		if(_index==2&&
			freArr[2]/2>190&&freArr[2]/2<204){ // ->197, 130~140format

			return freArr[2]/2;
		}

		if(_index==2&&
			dbArr[1]>dbArr[0]&&
			dbArr[2]-dbArr[1]>18&&
			freArr[2]>190&&freArr[2]<204){ // ->197, x124

			return freArr[3]/2;
		}

		if(_index==1&&
			dbArr[1]-dbArr[0]>18&&
			freArr[2]/2>240&&freArr[2]/2<255){ // ->247, x12

			return freArr[2]/2;
		}

		if(freArr[0]>130&&freArr[0]<160){ // ->80

			return freArr[0]/2;
		}

		if(_index==2&&
			vk3==3&&
			freArr[2]/2>100&&freArr[2]/2<120){ // ->110, 2369,x123

			return freArr[2]/2;
		}

		if(freArr[0]>65&&freArr[0]<75){ // ->71
			return freArr[0];
		}

		if(dbArr[0]-dbArr[1]>24&&
			dbArr[2]>dbArr[1]&&
			freArr[0]>220){ // ->247, 1x3

			return freArr[0];
		}

		if(dbArr[0]-dbArr[1]>12&&
			freArr[0]>220&&
			light>0.98){ // ->247, 1x3

			return freArr[0];
		}

		if(_index==1&&
			dbArr[1]-dbArr[0]>8&&
			freArr[1]>300&&
			heightArr[0]<15&&
			light>0.98&&
			refLength<6){ // ->329, x12

			return freArr[1];
		}

		if(_index==0){
			fre=freArr[0]/uk1;
			return fre;
		}
	}

	// 3469/369-10/369-11 ->330, valid ->110, valid ->246
	if(uk1==3&&uk2==4&&uk3==6&&
		vk3==9){

		if(refLength>9&&freArr[0]>800){ // ->330
			fre=freArr[0]/uk1;
		}
		else if(freArr[0]>100&&freArr[0]<120){ // ->110
			fre=freArr[2]/2;
		}
		else if(index1==0&&
			dbArr[2]>dbArr[1]&&
			freArr[0]>240&&freArr[0]<255){ // ->246, ❌82

			fre=freArr[2]/2;
		}
		
		return fre;
	}

	// 69-11+458 ->196
	if(uk1==6&&uk2==9&&
		vk1==4&&vk2==5&&vk3==8&&
		index1==1&&
		freArr[1]>190&&freArr[1]<205){

		return freArr[3]/2;
	}

	// 3467/3468 ->82, 70~90
	if(uk1==3&&uk2==4&&uk3==6&&
		(vk3==7||vk3==4)){

		if(freArr[0]>210&&freArr[0]<270){

			return freArr[0]/uk1;
		}
	}

	// 346 ->246, ❌82,  ->67, ->196
	if(uk1==3&&uk2==4&&uk3==6){
		if(index1==0&&
			dbArr[2]>dbArr[1]&&
			freArr[0]>240&&freArr[0]<255){ // ->246, ❌82

			return freArr[2]/2;
		}

		if(index1==0&&
			dbArr[1]>dbArr[2]&&
			freArr[0]>195&&freArr[0]<225){ // ->67, 65~75

			return freArr[0]/uk1;
		}

		if(index1==2&&
			freArr[2]>190&&freArr[2]<205&&
			vk3!=7){

			return freArr[2];
		}
	}

	// 679-12 ->valid 110, 2x34
	if(uk1==6&&uk2==7&&uk3==9&&
		vk3==12){

		if(index1==0&&
			freArr[0]>200&&freArr[0]<240){

			return freArr[0]/2;
		}
	}

	// 3679 ->valid-110
	if(uk1==3&&uk2==6&&uk3==7&&
		vk3==9){

		if(dbArr[0]>dbArr[1]&&
			dbArr[1]>dbArr[2]&&
			(dbArr[3]>dbArr[2]||
				dbArr[1]-dbArr[2]>12)){

			fre=freArr[0];
			return fre;
		}
	}

	// 367 ->valid-210, ->valid 197, ->110, ->247, 12x ->329, 12x
	if(uk1==3&&uk2==6&&uk3==7){
		if(dbArr[0]-dbArr[2]>18&&
			freArr[0]>190&&freArr[0]<204){ // valid ->197

			return freArr[0];
		}

		if(dbArr[1]-dbArr[2]>18&&
			freArr[1]>200&&freArr[1]<240){ // ->110

			return freArr[1]/2;
		}

		if(dbArr[0]>dbArr[1]&&
			dbArr[1]>dbArr[2]&&
			freArr[0]>100&&freArr[0]<120&&
			refLength>6){ // ->110 >6, 36???

			return freArr[1]/2;
		}

		if(dbArr[0]-dbArr[1]>18&&
			dbArr[1]>dbArr[2]&&
			dbArr[1]>dbArr[3]&&
			freArr[0]>220&&
			refLength<5){ // ->247, 12x

			return freArr[0];
		}

		if(dbArr[0]-dbArr[2]>18&&
			dbArr[1]-dbArr[2]>10&&
			dbArr[1]-dbArr[3]>10&&
			freArr[0]>220&&
			light>0.98&&
			refLength<5){ // ->247, 12x

			return freArr[0];
		}

		if(dbArr[0]-dbArr[2]>20&&
			dbArr[1]-dbArr[2]>18&&
			freArr[0]>300&&
			light>0.98&&
			refLength<6){ // ->329, 12x 

			return freArr[0];
		}

		if(dbArr[0]>dbArr[1]&&
			dbArr[1]-dbArr[2]>20){ // ->246, 12x

			fre=freArr[1]/2;
			return fre;
		}

		if(dbArr[0]>dbArr[1]&&
			dbArr[1]>dbArr[2]&&
			dbArr[1]>dbArr[3]){

			fre=freArr[0]/uk1;
			return fre;
		}
	}

	// 367 ->low-40, valid ->197, ->110
	if(uk1==3&&uk2==6&&uk3==7){
		if(dbArr[0]-dbArr[2]>18&&
			freArr[0]>190&&freArr[0]<204){ // valid ->197

			return freArr[0];
		}

		if(refLength<5&&
			dbArr[0]>dbArr[1]&&dbArr[0]>dbArr[2]){

			fre=freArr[0]/uk1;
			return fre;
		}

		if((vk3==9||(freArr[3]>300&&freArr[3]<360))&&
			index1==1&&
			freArr[1]>200&&freArr[1]<240){ // valid ->110

			return freArr[0];
		}
	}
	else if(uk1==3&&uk2==6&&uk3==8){ // 368 ->low-40, ->110, 3689,12x3
		if(dbArr[0]-dbArr[2]>18&&
			freArr[0]>190&&freArr[0]<204){ // valid ->197

			return freArr[0];
		}

		if(index1==1&&
			dbArr[1]-dbArr[2]>18&&
			freArr[1]/2>100&&freArr[1]/2<120){ // ->110, 3689,12x3

			return freArr[1]/2;
		}

		if(refLength<5&&
			dbArr[0]>dbArr[1]&&
			dbArr[0]>dbArr[2]&&
			dbArr[1]-dbArr[2]<12){

			fre=freArr[0]/uk1;
			return fre;
		}
	}

	// 356 ->330
	if(uk1==3&&uk2==5&&uk3==6){
		if(freArr[0]/3>300&&freArr[0]/3<360&&
			refLength>5){

			return freArr[0]/3;
		}
	}

	// valid 1246 ->nex, valid ->110, valid ->74, valid ->80, valid ->246
	if(uk1==1&&uk2==2&&uk3==4&&
		vk1==1&&vk2==2&&vk3==3){
		
		if(dbArr[1]>dbArr[2]&&
			dbArr[2]>dbArr[3]&&
			dbArr[3]>dbArr[0]&&
			freArr[0]>50&&freArr[0]<60){ // valid ->110

			return freArr[2]/2;
		}
		else if(freArr[0]<80){ // area ->74, valid ->110
			fre=freArr[0];

			if(freArr[0]<60){ // freArr[0]>50&&
				fre=freArr[2]/2;
			}

			return fre;
		}
		else{
			if(dbArr[0]-dbArr[1]>15){ // valid-110
				return freArr[0];
			}
		}

		if(index1==1&&
			freArr[0]>75&&freArr[0]<90){

			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[index1]/2, freArr[index1],1);
			if(flag){
				return freArr[index1]/2;
			}
		}

		if(index1==1&&
			freArr[1]>240&&freArr[1]<255){ // valid ->246

			return freArr[2]/2;
		}

		if(valid){
			*valid=1;
		}
		
		return 0;
	}

	// 1369
	if(uk1==1&&uk2==3&&uk3==6&&
		vk1==1&&vk2==2&&vk3==3){

		int _index=0;

		_index=__arr_maxIndex(dbArr+1, length-1);
		if(dbArr[0]-dbArr[_index+1]>10){
			return freArr[0];
		}
	}

	// valid 234 ->cut_valid low50~60, ->80
	/***
		1. 468 ->58 ->116
		2. 234 ->110 ->55
		3. 234 ->146 ->73
	****/
	if(uk1==2&&uk2==3&&uk3==4){
		if(freArr[1]<180&&freArr[1]>150){

			/***
				dbArr[0]>dbArr[2]&&
				dbArr[0]>dbArr[1]&&
				dbArr[2]>dbArr[3]&&
				dbArr[1]>dbArr[3]&&
			****/	

			if((fabsf(dbArr[0]-dbArr[1])<10||
				fabsf(dbArr[2]-dbArr[1])<10)&&
				fabsf(dbArr[0]-dbArr[2])<15&&
				((dbArr[1]-dbArr[3]>2&&heightArr[1]>15)||
					fabsf(freArr2[0]*2-freArr[0])<5||
					(dbArr[1]>dbArr[3]&&
						fabsf(dbArr[0]-dbArr[1])<6))){

				return freArr[0]/2;
			}

			if((fabsf(dbArr[0]-dbArr[1])<10||
				fabsf(dbArr[2]-dbArr[1])<10)&&
				fabsf(dbArr[0]-dbArr[2])<15&&
					dbArr[0]>dbArr[1]&&
					dbArr[2]>dbArr[1]&&
					dbArr[1]-dbArr[3]>3){

				return freArr[0]/2;
			}

			// if(fabsf(freArr2[0]*2-freArr[0])<5){
			// 	fre=freArr2[0];
			// }
			// else
			{
				int _flag=0;

				_flag=__queue_query(freArr3,dbArr3,heightArr3, refLength,freArr[0]/2);
				if(_flag){
					fre=freArr[0]/2;
				}
				else{
					fre=freArr[0];
				}

				// fre=__queue_cutValid(freArr,dbArr,length,0,1,
				// 					freArr2,dbArr2,length2,
				// 					freArr3,dbArr3,refLength);
			}
			
			return fre;
		}
		else if(freArr[0]>200&&freArr[0]<240){ // valid-110
			fre=freArr[0]/uk1;
			return fre;
		}
		// else if(freArr[0]<180&&freArr[0]>120){ // valid-74
		// 	fre=freArr[0]/uk1;
		// 	return fre;
		// }

		if(index1==2&&
			dbArr[0]>dbArr[1]&&
			freArr[0]>75&&freArr[0]<90){ // ->80,2x4

			return freArr[2]/2;
		}

		if((index1==0||index1==1)&&
			dbArr[index1]>dbArr[2]&&
			dbArr[index1]>dbArr[3]&&
			freArr[0]>150&&freArr[0]<180){ // ->80, 234

			return freArr[0]/2;
		}

	}
	else if(vk1==2&&vk2==3){ // valid ->80, ->197,x236(146-12),max-3
		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if(_index==1&&
			dbArr[0]-dbArr[2]<3&&
			freArr[1]>120&&freArr[1]<180){

			return freArr[1]/2;
		}

		if(uk2==4&&
			dbArr[1]-dbArr[0]>18&&
			freArr[1]>120&&freArr[1]<180){

			return freArr[1]/2;
		}

		if((_index==1||_index==2)&&
			freArr[1]/2>190&&freArr[1]/2<204){ // ->197

			return freArr[1]/2;
		}
	}
	
	//  valid ->71, valid ->110, valid ->197, 
	if(index1==1){
		int ts1=0,ts2=0;
		int tk1=0,tk2=0,tk3=0;

		int k1=0;
		int k2=0;

		int _index=0;

		_index=__arr_maxIndex(dbArr, length);

		__queue_fre3(freArr[1],freArr[2],freArr[3],
					&ts1,&ts2,
					&tk1,&tk2,&tk3);

		__queue_fre2(freArr[1],freArr[2],
					&k1,&k2);

		if((tk1==1||k1==1)&&
			dbArr[1]>dbArr[0]&&
			freArr[1]>190&&freArr[1]<204){ // x-1nn ->valid-197

			if(fabsf(freArr[0]*2-freArr[1])<4){
				fre=__queue_cutValid(freArr+1,dbArr+1,length-1,0,1,
									freArr2,dbArr2,length2,
									freArr3,dbArr3,refLength);

				return fre;
			}
			else{
				if(dbArr[1]>dbArr[2]&&
					dbArr[1]>dbArr[3]){

					return freArr[1];
				}
			}
		}

		if(tk1==2&&tk2==3&&tk3==4){ // x-234 ->valid-197/2
			if(freArr[1]>190&&freArr[1]<204&&
				_index==1&&
				dbArr[3]-dbArr[2]<6&&
				heightArr[2]>18){

				return freArr[1]/2;
			}
		}

		if(tk1==2&&tk2==3&&tk3==4){ // x-234 ->valid-110
			if(freArr[1]>210&&freArr[1]<230){
				return freArr[1]/2;
			}
		}

		// valid ->71
		if(k1==3&&k2==4&&
			freArr[1]>195&&freArr[1]<225&&
			fabsf(freArr[1]/3*4-freArr[2])<4){ 

			if(index1==1&&
				dbArr[1]-dbArr[2]>24&&
				freArr[0]>95&&freArr[0]<103){

				return freArr[1];
			}

			return freArr[1]/k1;
		}
		
		if(tk1==3&&tk2==4&&tk3==6){ // x-346 ->cut_valid-58
			if(dbArr[2]>dbArr[0]&&
				dbArr[2]>dbArr[1]&&
				dbArr[3]>dbArr[0]&&
				dbArr[3]>dbArr[1]&&
				freArr[1]<180&&freArr[1]>150){

				fre=freArr[1]/tk1;
				return fre;
			}
		}
	}

	// valid 2346 ->nex & area-74, valid ->110, valid ->197, valid ->80, valid ->147
	if(uk1==2&&uk2==3&&uk3==4&&
		vk1==3&&vk2==4&&vk3==6){

		int _index=0;

		_index=__arr_maxIndex(dbArr, length);

		if(refLength>6){ // ->147, 2346,x3x6-9; 2346,2x46
			int _index1=0;

			if(_index==3&&
				freArr[2]>280&&freArr[2]<310&&
				dbArr[0]-dbArr[1]>12&&
				dbArr[2]-dbArr[1]>12){ // 2x46 ->1x23

				return freArr[2]/2;
			}

			_index1=__arr_maxIndex(dbArr2+4, length2-4)+4;
			if(freArr[3]>280&&freArr[3]<310&&
				freArr2[3]>280&&freArr2[3]<310&&
				dbArr[3]>dbArr[2]){ // x3x6 ->x1x2-3

				float _fre=0;
				int k1=0,k2=0;

				if(_index==1){
					return freArr[3]/2;
				}

				if(freArr2[4]>420&&freArr2[4]<465){

					_fre=freArr2[4];
				}
				else if(freArr2[5]>420&&freArr2[5]<465){
					
					_fre=freArr2[5];
				}
				
				if(_fre){
					__queue_fre2(freArr[3],_fre,
								&k1,&k2);

					if(k1==2&&k2==3){
						return freArr[3]/2;
					}
				}
			}
		}

		if(dbArr[1]-dbArr[0]>12&&
			dbArr[1]-dbArr[2]>12&&
			dbArr[3]-dbArr[0]>12&&
			dbArr[3]-dbArr[2]>12){ // valid ->196, 2346,x3x6

			return freArr[3]/2;
		}

		if(_index==3||
			(_index==0&&dbArr[0]-dbArr[3]<2)){ 

			fre=freArr[0]/uk1;
			if(fre>60&&fre<80){
				return fre;
			}
		}
		else{
			if(_index<=1&&
				fabs(dbArr[0]-dbArr[1])<4){

				fre=freArr[0]/2;
				return fre;
			}
		}

		if(dbArr[0]>dbArr[2]&&
			dbArr[0]>dbArr[3]&&
			dbArr[1]>dbArr[2]&&
			dbArr[1]>dbArr[3]){

			if(fabsf(freArr[0]/2*3-freArr[1])<4&&
				freArr[0]>210&&freArr[0]<230){ // ->110
				
				return freArr[0]/2;
			}
		}

		if(!index1&&
			freArr[0]>150&&freArr[0]<170
			&&dbArr[1]>dbArr[2]){ // ->80

			return freArr[0]/2;
		}

		if(index1==2&&
			freArr[2]>150&&freArr[2]<170&&
			dbArr[0]>dbArr[1]&&
			dbArr[2]-dbArr[1]>15){ // ->80

			return freArr[2]/2;
		}

		if(index1==3&&
			freArr[2]>150&&freArr[2]<170&&
			dbArr[2]>dbArr[0]&&
			dbArr[0]>dbArr[1]){ // ->80

			return freArr[2]/2;
		}

		if(dbArr[0]-dbArr[1]>18&&
			freArr[0]>190&&freArr[0]<204){ // ->197

			return freArr[0];
		}

		if(index1==3&&
			freArr[index1]>230&&
			freArr[index1]<260&&
			dbArr[2]>dbArr[0]&&
			dbArr[2]>dbArr[1]){ // ->80

			return freArr[2]/2;
		}

		if(valid){
			*valid=1;
		}
		
		return 0;
	}

	// 1237 ->100
	if(uk1==1&&uk2==2&&uk3==3&&
		vk1==2&&vk2==3&&vk3==7){

		if(index1==1&&
			dbArr[0]>dbArr[2]&&
			dbArr[0]&&dbArr[3]){

			int _flag=0;

			for(int i=0;i<refLength;i++){
				if(fabsf(freArr[3]-freArr3[i])<2){
					_flag=1;
					break;
				}
			}

			if(_flag){
				return freArr[1]/2;
			}
		}
	}

	// valid 123/124 ->valid 58, ->valid 80, ->valid 197, ->valid 110, ->valid 147 ,1236,xx12
	if(uk1==1&&uk2==2&&(uk3==3||uk3==4)){
		if(index1==1&&
			freArr[0]>60&&freArr[0]<85){ // ->80

			return freArr[1]/2;
		}

		if(uk3==3&&
			freArr[0]>190&&freArr[0]<204&&
			dbArr[1]-dbArr[0]<3){ // ->197

			return freArr[0];
		}

		if(!index1&&uk3==3&&
			freArr[1]>200&&freArr[1]<240){ // ->110

			return freArr[1]/2;
		}

		if(freArr[index1]>150&&freArr[index1]<170&&
			refLength>5){ // ->80

			float _fre=0;

			int us1=0,us2=0;
			int uk1=0,uk2=0,uk3=0;

			int k1=0,k2=0;

			for(int i=0;i<refLength-2;i++){
				if(freArr3[i]>freArr[index1]){
					_fre=__queue_fre3(freArr3[i],freArr3[i+1],freArr3[i+2],
									&us1,&us2,
									&uk1,&uk2,&uk3);

					if(us1==1&&us2==1&&
						freArr[index1]>_fre){

						__queue_fre2(_fre,freArr[index1],
									&k1,&k2);

						if(k1==1&&k2==2){
							return freArr[index1]/2;
						}
					}
				}
			}
		}

		if(uk3==3&&
			index1==2&&
			freArr[2]>280&&freArr[2]<310){ // ->147

			int flag=0;

			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[2]/2, freArr[2],0);
			if(flag){
				return freArr[2]/2;
			}
		}

		if(freArr[1]>190&&freArr[1]<204&&
			dbArr[0]-dbArr[1]<6&&
			refLength>5){

			fre=__queue_cutValid(freArr+1,dbArr+1,length-1,0,1,
								freArr2,dbArr2,length2,
								freArr3,dbArr3,refLength);

			return fre;
		}

		if(freArr[0]>50&&freArr[0]<60&&
			dbArr[1]>dbArr[2]&&
			dbArr[2]-dbArr[0]>12){

			int _flag=0;

			_flag=__queue_query(freArr3,dbArr3,heightArr3, refLength,freArr[0]);
			if(_flag){
				fre=freArr[1]/2;
			}
			else{
				fre=freArr[2]/2;
			}
		}

		if(dbArr[0]>dbArr[1]&&
			dbArr[1]>dbArr[2]&&
				dbArr[1]>dbArr[3]&&
				refLength>5){

			int k1=0,k2=0;

			fre=__queue_cutValid(freArr,dbArr,length,0,0,
								freArr2,dbArr2,length2,
								freArr3,dbArr3,refLength);

			__queue_fre2(fre, freArr[0], &k1, &k2);
			if(k1==1&&k1==k2){
				fre=freArr[1]/2;
			}

			return fre;
		}
	}

	// area ->74 1346/146/346; 5string ->110, 346; 4string ->147, 346; ->196
	if(uk1==1&&uk2==3&&uk3==4&&
		vk3==6){

		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if(_index==3){
			fre=freArr[0]/uk1;
			if(freArr[3]>190&&freArr[3]<205){ // ->196
				return freArr[3];
			}

			return fre;
		}

		if(!_index&&
			freArr[1]/3>105&&freArr[1]/3<115){ // valid ->110

			return freArr[1]/3;
		}
	}
	else if((uk1==1&&uk2==4&&uk3==6)||
		(uk1==3&&uk2==4&&uk3==6)){

		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if(_index==2&&
			dbArr[1]>dbArr[0]&&
			dbArr[1]>dbArr[3]){

			if(freArr[2]>190&&freArr[2]<204&&
				dbArr[2]-dbArr[1]>15){ // ->196

				return freArr[2];
			}

			if(freArr[1]/2>190&&freArr[1]/2<204){ // ->196

				return freArr[1]/2;
			}

			fre=freArr[0]/uk1;
			if(uk1==3&&fre<65){ // 74-9
				fre=freArr[1]/2;
			}

			return fre;
		}
	}

	// area ->74 1457/2456, ->110, 2456
	if(uk1==1&&uk2==4&&uk3==5&&
		vk3==7){

		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if(_index==2){
			return freArr[0];
		}
	}
	else if(uk1==2&&uk2==4&&uk3==5&&
		vk3==6){

		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if((_index==3||
			(_index==0&&dbArr[0]-dbArr[3]<2))&&
			freArr[0]>120&&freArr[0]<160){

			return freArr[0]/2;
		}

		if(freArr[0]>200&&freArr[0]<240){
			return freArr[0]/2;
		}
	}

	// valid ->110 ,1x23-4, 1x36, 12x3; ->80
	if(!uk1){
		if(dbArr[0]>dbArr[2]&&
			dbArr[2]>dbArr[1]&&
			dbArr[2]>dbArr[3]){ // 1x23-4

			_fre=__queue_fre3(freArr[0],freArr[2],freArr[3],
							&us1,&us2,
							&uk1,&uk2,&uk3);

			if(uk1==1){
				if(freArr[2]/uk2>105&&freArr[2]/uk2<115){
					return freArr[0]/uk1;
				}
			}

			if(freArr[2]>210&&freArr[2]<230){
				int k1=0;
				int k2=0;

				_fre=__queue_fre2(freArr[2],freArr[3],
								&k1,&k2);

				if(k1==2&&k2==3){
					return freArr[2]/2;
				}
			}
		}

		if(dbArr[0]>dbArr[1]&&
			dbArr[1]>dbArr[2]&&
			dbArr[1]>dbArr[3]){

			int k1=0;
			int k2=0;

			_fre=__queue_fre2(freArr[0],freArr[2],
							&k1,&k2);

			if(k1==1){
				if(freArr[2]/uk2>105&&freArr[2]/uk2<115){
					return _fre; 
				}
			}
		}

		if(dbArr[0]>dbArr[2]&&
			dbArr[1]>dbArr[2]&&
			dbArr[3]>dbArr[2]){ // 12x3

			_fre=__queue_fre3(freArr[0],freArr[1],freArr[3],
							&us1,&us2,
							&uk1,&uk2,&uk3);

			if(uk1==1&&
				freArr[1]>210&&freArr[1]<230){

				return freArr[0]/uk1;
			}
		}
	}
	else{
		if(dbArr[0]>dbArr[2]&&
			dbArr[0]>dbArr[3]&&
			dbArr[1]>dbArr[2]&&
			dbArr[1]>dbArr[3]){ // 12xx

			if(uk1==1&&uk2==2&&
				freArr[1]>210&&freArr[1]<230){

				return freArr[1]/2;
			}
		}

		if(dbArr[0]>dbArr[2]&&
			dbArr[1]>dbArr[2]&&
			dbArr[3]>dbArr[2]){ // 12x3

			_fre=__queue_fre3(freArr[0],freArr[1],freArr[3],
							&us1,&us2,
							&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==2&&
				freArr[1]>210&&freArr[1]<230){

				return freArr[0];
			}
		}
	}

	// n124, ->330
	if(vk1==1&&vk2==2&&vk3==4&&
		index1==3){

		if(freArr[1]>220&&freArr[1]<360){
			return freArr[2]/2;
		}
	}

	// valid ->80,{77,90} 160-max, 12/23/query2
	if(freArr[index1]>154&&
		freArr[index1]<180&&
		refLength>3){

		int flag=0;

		int k1=0;
		int k2=0;

		int _cutLen=0;

		// 256-8
		if(index1==0&&
			uk1==2&&uk2==5&&uk3==6){

			return freArr[index1]/2;
		}

		// 12
		// if(index1>0){
		// 	for(int i=0;i<index1;i++){
		// 		__queue_fre2(freArr[i],freArr[index1],
		// 					&k1,&k2);

		// 		if(k1==1&&k2==2&&
		// 			fabsf(freArr[i]*2-freArr[index1])<5){

		// 			return freArr[index1]/2;
		// 		}
		// 	}
		// }

		// 23 ->query2 ???
		if(index1<2){
			__queue_fre2(freArr[index1],freArr[index1+1],
						&k1,&k2);

			if(k1==2&&k2==3){
				// if(dbArr[index1]-dbArr[index1+1]<15){
				// 	flag=1;
				// }

				flag=1;
				if(flag){
					return freArr[index1]/2;
				}
			}
			else{
				if(fabsf(freArr[index1]/2-freArr[index1+1]/3)<5){  // error<15
					// if(dbArr[index1]-dbArr[index1+1]<15){
					// 	flag=1;
					// }

					flag=1;
					if(flag){
						return freArr[index1]/2;
					}
				}
			}
		}

		// query2
		flag=1;
		if(index1==1&&
			dbArr[index1]-dbArr[index1+1]>20){

			__queue_fre2(freArr[index1],freArr[index1+2],
						&k1,&k2);

			if(k1==1&&k2==2){
				flag=0;
			}
		}

		if(flag){
			_cutLen=__arr_cut(freArr3, refLength, freArr[index1]*4+10);
			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[index1]/2, freArr[index1],1);
			if(flag){
				return freArr[index1]/2;
			}
		}
	}

	// valid ->80,{75,90}, 27-11
	if(freArr[0]>150&&
		freArr[0]<180&&
		refLength>3){

		int flag=0;

		__queue_fre3(freArr[0],freArr[1],freArr[2],
					&us1,&us2,
					&uk1,&uk2,&uk3);

		if(uk1==2&&uk2==3){
			return freArr[0]/2;
		}
		else if(uk1==1){
			flag=__queue_query2(freArr3, dbArr3, heightArr3, refLength, 0, freArr[0]/2, freArr[0],1);
			if(flag){
				return freArr[0]/2;
			}
		}

		if(dbArr[0]>dbArr[2]&&
			dbArr[0]>dbArr[3]&&
			fabsf(freArr[0]/2*7-freArr[1])<4){ // 27-11

			return freArr[0]/2;
		}
	}

	// valid ->329, 1236
	if(uk1==1&&uk2==2&&uk3==3&&
		vk3==6&&
		index1==2&&
		freArr[index1]>315&&freArr[index1]<345){

		int flag=0;
		int c1=0;

		flag=__queue_four(freArr3, dbArr3, heightArr3, refLength, freArr[index1]);
		if(flag){
			c1=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr[index1]*4+20, freArr[index1], 1);
			if(c1){
				return freArr[index1];
			}
		}
	}

	// valid ->110, ->147, 12/1; ->6string ->80
	{
		float _fre=0;

		int k1=0;
		int k2=0;

		// dB desc
		__vcorrsort1(dbArr, freArr,heightArr,indexArr,length, 1);
		// fre asc 2
		__vcorrsort1(freArr, dbArr,heightArr,indexArr,2, 0);

		_fre=__queue_fre2(freArr[0],freArr[1],
						&k1,&k2);

		if(k1==2&&k2==3&&
			fabsf(freArr[0]/2*3-freArr[1])<4){

			if(freArr[0]>210&&freArr[0]<230){ // ->110, ->329
				
				if(dbArr[1]-dbArr[0]>6&&
					light>0.98&&
					heightArr[0]<15&&
					refLength<6){ // ->329, x14

					return freArr[1];
				}
				else if(dbArr[1]-dbArr[0]>12&&
					heightArr[0]<8&&
					refLength<=3){ // ->329, 

					return freArr[1];
				}
				else{
					return _fre;
				}
			}
			else if(freArr[0]>140&&freArr[0]<180&&
					fabsf(dbArr[0]-dbArr[1])<12){ // ->80 

				return _fre;
			}
		}

		if(k1==1&&k2==2&&
			fabsf(freArr[0]*2-freArr[1])<4){

			if(freArr[0]>130&&freArr[0]<160){ // ->146, 80 ???
				return freArr[1]/2;
			}
			else if(freArr[0]>60&&freArr[0]<85){ // ->80
				return freArr[1]/2;
			}
			else if(freArr[0]>190&&freArr[0]<204){ // ->197
				return freArr[0];
			}
		}
	}

	// dB desc
	__vcorrsort1(dbArr, freArr,heightArr,indexArr,length, 1);
	// fre asc 3
	__vcorrsort1(freArr, dbArr,heightArr,indexArr,3, 0);

	_fre=__queue_fre3(freArr[0],freArr[1],freArr[2],
					&us1,&us2,
					&uk1,&uk2,&uk3);

	if(!uk1){ // valid ->110, 1x3 (x=2)
		float _fre=0;

		int k1=0;
		int k2=0;

		_fre=__queue_fre2(freArr[0],freArr[2],
						&k1,&k2);

		if(k1==1&&k2==3&&
			fabsf(freArr[0]*3-freArr[2])<4&&
			freArr[0]>100&&freArr[0]<200){

			if(fabsf(freArr[0]*2-freArr[1])<10){
				return _fre;
			}
		}
	}

	if(uk1==1&&uk2==3&&uk3==4){
		if(freArr[0]>100&&freArr[0]<120){

			return freArr[1]/3;
		}
	}
				
	if(us1==1&&us1==us2){
		if(fabsf(_fre*uk2-freArr[1])<5&&
			fabsf(_fre*uk3-freArr[2])<5){

			fre=_fre;

			// valid 1:x:2
			index1=__arr_maxIndex(dbArr, 3);
			if(index1==0){
				if(uk1==2&&2*uk1==uk3&&
					length2>=4){ // valid-110 ->234-5, 234-7, ->247,234->1x2

					__queue_fre3(freArr2[1],freArr2[2],freArr2[3],
								&vs1,&vs2,
								&vk1,&vk2,&vk3);

					if(vk1==3&&vk2==4&&(vk3==5||vk3==7)){
						return fre;
					}
				}

				if(dbArr[0]-dbArr[1]>20&&
					dbArr[2]-dbArr[1]>10&&
					freArr[0]>220){ // ->247,234->1x2

					return freArr[2]/2;
				}
				
				if(2*uk1==uk3&&refLength>5){
					if(valid){
						*valid=3;
					}

					return 0;
				}
			}

			// valid 2:x:3
			if(uk1==4&&uk3==6&&refLength>5){
				if(valid){
					*valid=3;
				}

				return 0;
			}
		}
	}
	else{
		// 236 ->valid-210 low65~75
		if(uk1==2&&uk2==3&&uk3==6){
			if(dbArr[1]>dbArr[2]&&
				dbArr[2]>dbArr[0]&&
				dbArr[0]>dbArr[3]&&
				freArr[0]<150&&freArr[0]>130){

				fre=freArr[0]/uk1;
				return fre;
			}
		}

		// 236 ->valid-210
		if(uk1==2&&uk2==3&&uk3==6&&refLength>5){
			int _index=0;

			_index=__arr_maxIndex(dbArr, length);
			if(_index==1){
				if(dbArr[1]>dbArr[2]&&
					dbArr[2]>dbArr[0]){

					if(valid){
						*valid=1;
					}
					
					return 0;
				}
			}
		}
	}

	// ->71, aug format 3times, 1-236,34, valid ->196, 185~205
	if(!fre&&refLength<5){
		if(vk1==2&&vk2==3&&vk3==6){
			int _index=0;

			_index=__arr_maxIndex(dbArr, length);
			if(_index==2){
				if(dbArr[2]-dbArr[1]>15&&
					dbArr[1]>dbArr[0]){

					fre=freArr[0];
				}
			}
		}
		else if(!vk1&&!uk1){
			int _index=0;

			_index=__arr_maxIndex(dbArr, length);

			if(_index==1&&
				dbArr[1]-dbArr[0]>12&&
				dbArr[1]-dbArr[2]>20&&
				dbArr[1]-dbArr[3]>20){

				int k1=0;
				int k2=0;

				__queue_fre2(freArr[0], freArr[1], &k1, &k2);
				if(fabsf(freArr[0]/k1*k2-freArr[1])>5&&
					freArr[1]>190&&freArr[1]<200){
					
					return freArr[1];
				}
			}

			if(_index==1){
				if(dbArr[1]-dbArr[0]>15&&
					dbArr[1]-dbArr[2]>15&&
					dbArr[1]-dbArr[3]>15){

					int k1=0;
					int k2=0;

					__queue_fre2(freArr[1], freArr[2], &k1, &k2);
					if(k1==3&&k2==4){
						if(dbArr[1]-dbArr[2]>30){ // ->x1xx
							return freArr[1];
						}

						return freArr[1]/3;
					}
				}
			}
		}
	}

	// boundary ->string1/3/4/5/6
	if(!fre){
		float _fre=0;

		int _index=0;
		int k1=0,k2=0;
		int tk1=0,tk2=0;

		// dB desc
		__vcorrsort1(dbArr, freArr,heightArr,indexArr,length, 1);
		// fre asc 3
		__vcorrsort1(freArr, dbArr,heightArr,indexArr,length, 0);

		_index=__arr_maxIndex(dbArr, length);

		__queue_fre3(freArr[0],freArr[1],freArr[2],
					&us1,&us2,
					&uk1,&uk2,&uk3);

		// ->1string 330
		if(_index&&freArr[_index]>520){
			_fre=__queue_fre2(freArr[_index-1],freArr[_index],
							&k1,&k2);

			if(_index>=2&&
				freArr[_index-2]>140&&freArr[_index-2]<155){ // 12n,13n
				
				int flag=0;

				__queue_fre2(freArr[_index-2],freArr[_index-1],
							&tk1,&tk2);

				if(tk1==1&&(tk2==2||tk2==3)){
					flag=__queue_query2(freArr3,dbArr3,heightArr3,refLength,0,freArr[index1-2],freArr[index1-2]*2,0);
					if(flag){
						return freArr[index1-1]/tk2;
					}
				}
			}
			else if(_fre>280&&_fre<310){
				int flag=0;

				flag=__queue_query2(freArr3,dbArr3,heightArr3,refLength,0,_fre/2,_fre,0);
				if(flag){
					return _fre/2;
				}
			}

			if(k1==1){
				fre=_fre;
				return fre;
			}
		}

		// ->2string 247
		if(index1==3&&
			freArr[3]/2>240&&freArr[3]/2<255){

			__queue_fre2(freArr[index1-1],freArr[index1],
							&tk1,&tk2);

			if(tk1==1&&tk2==2&&
				fabsf(freArr[index1-1]*tk2-freArr[index1])<5){

				return freArr[index1]/2;
			}
		}

		// ->3string 197
		if(_index<3){
			_fre=__queue_fre2(freArr[_index],freArr[_index+1],
							&k1,&k2);
			
			if(k1==1&&
				freArr[_index]>190&&freArr[_index]<204){

				fre=freArr[_index];
				return fre;
			}
			else if(_index==1&&
				freArr[_index]>190&&freArr[_index]<204&&
				dbArr[1]-dbArr[2]>20&&
				dbArr[3]>dbArr[2]){
				
				_fre=__queue_fre2(freArr[1],freArr[3],
								&k1,&k2);
				if(k1==1){
					if(k2<5){
						return freArr[3]/k2;
					}
					else{
						return freArr[1];
					}
				}
			}
		}

		// ->6string 75~95
		_fre=__queue_fre2(freArr[0],freArr[1],
						&k1,&k2);
		
		if(_index==1){
			if(k1==1&&k2==2&&
				fabsf(freArr[0]*2-freArr[1])<5&&
				freArr[0]>70&&freArr[0]<90&&
				dbArr[0]-dbArr[2]<3&&
				dbArr[0]-dbArr[3]<3){

				fre=_fre;
				return fre;
			}
		}

		// ->4string 147 ->49/37/29 346/458/56-10, 1, x1x2
		if(_index==1){
			if(dbArr[1]-dbArr[0]>15){
				int us1=0,us2=0;
				int uk1=0,uk2=0,uk3=0;

				__queue_fre3(freArr[1],freArr[2],freArr[3],
							&us1,&us2,
							&uk1,&uk2,&uk3);
				
				if((uk1==3&&uk2==4&&uk3==6)||
					(uk1==4&&uk2==5&&uk3==8)||
					(uk1==5&&uk2==6&&uk3==10)){

					fre=freArr[1];
					return fre;
				}
			}

			if(freArr[1]>130&&freArr[1]<160){

				if(dbArr[1]-dbArr[2]>20&&
					dbArr[1]-dbArr[3]>20){

					return freArr[1];
				}
				else{ // 80 ???
					float _fre=0;

					int k1=0;
					int k2=0;

					_fre=__queue_fre2(freArr[1],freArr[3],
									&k1,&k2);
					if(k1==1){
						return freArr[1];
					}
				}
			}
		}

		// ->5string 110->37
		if(!_index&&
			uk1==3&&(uk2==4||uk2==5)&&uk3==6){

			if(dbArr[0]>dbArr[1]&&
				dbArr[2]>dbArr[1]){

				int us1=0,us2=0;
				int uk1=0,uk2=0,uk3=0;

				_fre=__queue_fre3(freArr[0],freArr[2],freArr[3],
								&us1,&us2,
								&uk1,&uk2,&uk3);

				if(uk1==1&&
					freArr[2]>210&&freArr[2]<230){

					fre=_fre;
					return fre;
				}
			}
		}
		else{
			if(_index){ // 23x
				if(dbArr[0]>dbArr[2]&&
					dbArr[0]>dbArr[3]&&
					dbArr[1]>dbArr[2]&&
					dbArr[1]>dbArr[3]){

					_fre=__queue_fre2(freArr[_index-1],freArr[_index],
									&k1,&k2);

					if(k1==2&&k2==3&&
						freArr[_index-1]>210&&freArr[_index-1]<230){

						if(index1==1&&
							refLength<=3&&
							heightArr[0]<8&&
							dbArr[1]-dbArr[0]>12){

							return 0;
						}

						fre=_fre;
						return fre;
					}
				}
			}
		}

		// ->6string 75~95
		_fre=__queue_fre2(freArr[0],freArr[1],
						&k1,&k2);
		
		if(k1==1&&k2==2&&
			fabsf(freArr[0]*2-freArr[1])<5&&
			freArr[0]<95&&
			dbArr[1]-dbArr[0]<12&&
			dbArr[0]>dbArr[2]&&
			dbArr[0]>dbArr[3]){

			fre=_fre;
			return fre;
		}
		else if(index1&&index1<3&&
			freArr[index1]>150&&freArr[index1]<170){ // 75~85

			_fre=__queue_fre2(freArr[index1],freArr[index1+1],
							&k1,&k2);

			if(k1==2&&k2==3){
				return freArr[index1]/2;
			}

			if(index1==1&&
				dbArr[0]>dbArr[2]&&
				dbArr[0]>dbArr[3]&&
				fabsf(freArr[0]-freArr[1]/2)<5){

				return freArr[1]/2;
			}

			if(index1==1&&
				dbArr[0]>dbArr[2]&&
				dbArr[0]>dbArr[3]){

				for(int i=0;i<refLength-2;i++){
					if(freArr3[i]>freArr[index1]){
						_fre=__queue_fre3(freArr3[i],freArr3[i+1],freArr3[i+2],
										&us1,&us2,
										&uk1,&uk2,&uk3);

						if(us1==1&&us2==1&&
							freArr[index1]>_fre){

							__queue_fre2(_fre,freArr[1],
										&k1,&k2);

							if(k1==1&&k2==2){
								return freArr[1]/2;
							}
						}
					}
				}
			}
		}
	}

	// ->valid 110, 12, 1<2; 1x23, 2>1>3
	if(!fre&&
		freArr[0]>105&&freArr[0]<115){

		int us1=0,us2=0;
		int uk1=0,uk2=0,uk3=0;

		int k1=0;
		int k2=0;

		// 1x24
		if(index1==2&&
			dbArr[0]>dbArr[1]&&
			dbArr[3]>dbArr[1]){

			__queue_fre3(freArr[0],freArr[2],freArr[3],
						&us1,&us2,
						&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==2&&
				fabsf(freArr[0]*2-freArr[2])<5){

				return freArr[2]/2;
			}
		}

		// 12nn
		if(index1==1&&
			dbArr[0]>dbArr[2]&&
			dbArr[0]>dbArr[3]){ // dbArr[1]-dbArr[0]>5&&

			__queue_fre2(freArr[0],freArr[1],
						&k1,&k2);

			if(k1==1&&k2==2&&
				fabsf(freArr[0]*2-freArr[1])<5){

				return freArr[1]/2;
			}
		}

		// 1x23
		if(dbArr[0]-dbArr[3]>12&&
			dbArr[2]-dbArr[3]>12){

			__queue_fre3(freArr[0],freArr[2],freArr[3],
						&us1,&us2,
						&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==2&&uk3==3&&
				fabsf(freArr[0]*2-freArr[2])<4&&
				fabsf(freArr[0]*3-freArr[3])<4){

				return freArr[2]/2;
			}
		}

		// 1x3
		if(dbArr[0]-dbArr[2]>12&&
			dbArr[2]-dbArr[3]>20){

			__queue_fre2(freArr[0],freArr[2],
						&k1,&k2);

			if(k1==1&&k2==3&&
				fabsf(freArr[0]-freArr[2]/3)<4){

				return freArr[2]/3;
			}
		}
	}

	// valid ->110, 220-max, 124/1x24
	if(!fre&&
		freArr[index1]>200&&
		freArr[index1]<240&&
		index1<3&&
		refLength>10){

		int flag=0;

		if(index1==0&&
			dbArr[0]-dbArr[1]>15){
			
			flag=1;
		}
		else{
			if(dbArr[index1]-dbArr[index1-1]>15&&
				dbArr[index1]-dbArr[index1+1]>15){

				flag=1;
			}
		}
		
		if(flag){
			flag=__queue_query2(freArr3,dbArr3,heightArr3,refLength,0,freArr[index1]/2,freArr[index1],0);
			if(flag){
				return freArr[index1]/2;
			}
		}

		if(index1==1&&
			fabsf(freArr[index1]/2-freArr[0])<5){ // 124

			flag=__queue_query2(freArr3,dbArr3,heightArr3,refLength,0,freArr[index1]/2,freArr[index1],1);
			if(flag){
				return freArr[index1]/2;
			}
		}

		if(index1==2&&
			dbArr[0]>dbArr[1]&&
			dbArr[2]>dbArr[1]&&
			fabsf(freArr[index1]/2-freArr[0])<5){ // 1x24

			flag=__queue_query2(freArr3,dbArr3,heightArr3,refLength,0,freArr[index1]/2,freArr[index1],1);
			if(flag){
				return freArr[index1]/2;
			}
		}
	}

	// valid ->110, 123, 2<1<3, 13-n, 1x23
	if(!fre&&
		((freArr[index1]>315&&freArr[index1]<345)||
			(freArr[index1]>105&&freArr[index1]<115))&&
		refLength>10){

		int us1=0,us2=0;
		int uk1=0,uk2=0,uk3=0;

		int k1=0;
		int k2=0;

		int flag=0;

		// 123
		__queue_fre3(freArr[0],freArr[1],freArr[2],
					&us1,&us2,
					&uk1,&uk2,&uk3);

		if(uk1==1&&uk2==2&&uk3==3){
			flag=__queue_query3(freArr3,dbArr3,heightArr3,refLength,0,freArr[2]/3,freArr[2],0);
			if(flag){
				return freArr[2]/3;
			}
		}

		// 13
		__queue_fre2(freArr[0],freArr[1],
					&k1,&k2);

		if(k1==1&&k2==3){
			flag=__queue_query3(freArr3,dbArr3,heightArr3,refLength,0,freArr[1]/3,freArr[1],0);
			if(flag){
				return freArr[1]/3;
			}
		}

		// 123
		__queue_fre3(freArr[0],freArr[2],freArr[3],
					&us1,&us2,
					&uk1,&uk2,&uk3);

		if(uk1==1&&uk2==2&&uk3==3){
			flag=__queue_query3(freArr3,dbArr3,heightArr3,refLength,0,freArr[3]/3,freArr[3],0);
			if(flag){
				return freArr[3]/3;
			}
		}
	}

	// valid 110, max-220, 2x35,12x5
	if(!fre&&
		freArr[index1]>200&&
		freArr[index1]<240&&
		refLength>5){

		int us1=0,us2=0;
		int uk1=0,uk2=0,uk3=0;

		int k1=0;
		int k2=0;

		if(index1==0){ // 2x35
			__queue_fre2(freArr[0],freArr[2],
						&k1,&k2);

			if(k1==2&&k2==3&&
				fabsf(freArr[0]/2*3-freArr[2])<5){

				return freArr[0]/2;
			}
		}

		if(index1==1){ // 12x5
			__queue_fre3(freArr[0],freArr[1],freArr[3],
						&us1,&us2,
						&uk1,&uk2,&uk3);


			if(uk1==1&&uk2==2&&
				fabsf(freArr[0]*2-freArr[1])<5&&
				fabsf(freArr[0]*uk3-freArr[3])<uk3*3){ // &&fabsf(freArr[0]*uk3-freArr[3])<uk3*3

				return freArr[1]/2;
			}
		}
	}

	// valid 80, 1x23,x236, 3>2>>1
	if(!fre&&
		freArr[index1]>230&&
		freArr[index1]<260&&
		index1>1){

		int k1=0;
		int k2=0;

		int flag=0;

		__queue_fre2(freArr[index1-1],freArr[index1],
					&k1,&k2);

		if((k1==2&&k2==3)||
			fabsf(freArr[index1-1]/2-freArr[index1]/3)<5){ // error 15

			if(index1==2&&
				dbArr[1]>dbArr[0]&&
				dbArr[1]>dbArr[3]){

				flag=1;
			}
			else if(dbArr[2]>dbArr[0]&&
				dbArr[2]>dbArr[1]){

				int _len;

				_len=refLength-1;
				for(int i=0;i<refLength;i++){
					if(freArr3[i]>1200){
						_len=i;
					}
				}

				if(_len>5){
					flag=1;
				}
			}

			if(flag){
				return freArr[index1-1]/2;
			}
		}
	}

	// valid 80, x23n
	if(!fre&&
		dbArr[1]>dbArr[0]&&
		dbArr[2]>dbArr[0]&&
		dbArr[1]>dbArr[3]&&
		dbArr[2]>dbArr[3]){

		int k1=0;
		int k2=0;

		__queue_fre2(freArr[1],freArr[2],
					&k1,&k2);

		if(k1==2&&k2==3&&
			freArr[1]>140&&freArr[1]<180&&
			fabsf(freArr[1]/2-freArr[2]/3)<3){

			return freArr[1]/2;
		}
	}

	// valid 147, x124, x245, 1x23
	if(!fre&&
		freArr[index1]>280&&
		freArr[index1]<310&&
		refLength>3){

		if(index1==2){
			__queue_fre3(freArr[index1-1],freArr[index1],freArr[index1+1],
						&us1,&us2,
						&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==2){
				return freArr[index1]/2;
			}

			__queue_fre3(freArr[0],freArr[index1],freArr[index1+1],
						&us1,&us2,
						&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==2&&uk3==3){ // 1x23
				return freArr[index1]/2;
			}
		}
		else if(index1==1){
			int _index=0;

			__queue_fre3(freArr[index1],freArr[index1+1],freArr[index1+2],
						&us1,&us2,
						&uk1,&uk2,&uk3);

			if(uk1==2&&uk2==4&&uk3==5){
				return freArr[index1]/2;
			}

			_index=__arr_maxIndex(dbArr2, length2);
			if(_index+2<length2){
				__queue_fre3(freArr2[_index],freArr2[_index+1],freArr2[_index+2],
						&us1,&us2,
						&uk1,&uk2,&uk3);

				if((uk1==2&&uk2==3)||
					(uk1==2&&uk2==4&&uk3==5)){

					return freArr2[_index]/2;
				}

			}
		}
	}

	// valid 196, x12n, 2-max
	if(!fre&&
		index1==2&&
		dbArr[1]>dbArr[0]&&
		dbArr[1]>dbArr[3]){

		int k1=0;
		int k2=0;

		__queue_fre2(freArr[1],freArr[2],
					&k1,&k2);

		if(k1==1&&k2==2&&
			fabsf(freArr[1]-freArr[2]/2)<8){

			int flag=0;

			flag=__queue_query2(freArr3,dbArr3,heightArr3,refLength,0,freArr[index1]/2,freArr[index1],1);
			if(flag){
				return freArr[index1]/2;
			}
		}
	}

	// valid 246, x123,3-max
	if(!fre&&
		index1==2&&
		freArr[2]/2>230&&
		refLength>12){

		__queue_fre3(freArr[1],freArr[2],freArr[3],
						&us1,&us2,
						&uk1,&uk2,&uk3);

		if(uk1==1&&uk2==2&&uk3==3){
			return freArr[2]/2;
		}
	}

	// valid->246, light=1, refLength>6, queue_multi
	if(!fre&&
		light>0.98&&
		refLength>6){

		float _fre1=0;
		float _fre2=0;

		int _num=2;
		int _subType=0;
		int _unionType=1;
		int _direction=0;

		int k1=0;
		int k2=0;

		_fre1=__queue_multi(freArr3,dbArr3,heightArr3,refLength,_num,_subType, _unionType, _direction);
		if(_fre1>230&&_fre1<255){
			// _fre2=freArr[index1];
			// if(_fre1>_fre2){
			// 	__queue_fre2(_fre2,_fre1,
			// 				&k1,&k2);

			// 	if(!k1){
			// 		fre=_fre1;
			// 	}
			// }

			fre=_fre1;
		}
		else if(_fre1>300&&_fre1<345){
			fre=_fre1;
		}

		if(!fre&&light>0.99){
			_unionType=2;
			_fre1=__queue_multi(freArr3,dbArr3,heightArr3,refLength,_num,_subType, _unionType, _direction);
			if(_fre1>300&&_fre1<345){
				fre=_fre1;
			}
		}

		if(!fre&&
			freArr[2]>240&&freArr[2]<255){

			int flag=0;
			int _index=0;

			flag=__queue_bear(freArr3, dbArr3, heightArr3, refLength, 1500, freArr[2], &_index);
			if(flag){
				return freArr[2];
			}
		}
	} 

	// valid 246,
	if(!fre&&
		refLength>9){

		float _fre1=0;

		int _num=2;
		int _subType=0;
		int _unionType=1;
		int _direction=0;

		_fre1=__queue_multi(freArr3,dbArr3,heightArr3,refLength,_num,_subType, _unionType, _direction);
		if(_fre1>230&&_fre1<255){
			fre=_fre1;
		}
	}

	// valid 329
	if(!fre&&
		freArr[index1]>300&&freArr[index1]<360){

		int flag=0;
		int c1=0;

		flag=__queue_four(freArr3, dbArr3, heightArr3, refLength, freArr[index1]);
		if(flag){
			c1=__queue_count(freArr3, dbArr3, heightArr3, refLength, 0, freArr[index1]*4+20, freArr[index1], 1);
			if(c1){
				fre=freArr[index1];
			}
		}
	}

	// valid 329
	if(!fre&&
		index1==1&&
		freArr[index1]>300&&freArr[index1]<360&&
		freArr[3]>2000&&
		refLength>4){

		int k1=0;
		int k2=0;

		__queue_fre2(freArr[1], freArr[2], &k1, &k2);
		if(k1==1&&k2==2){
			__queue_fre2(freArr3[3], freArr3[4], &k1, &k2);
			if(k1+1==k2&&
				fabsf(freArr3[3]/k1-freArr[1])<10){

				fre=freArr[2]/2;
			}
		}
	}

	return fre;
}

/***
	1. 3~4 ->1:1 strict max-3~4
	2. 1:2 -> strict max-2 ,use class 3 ???
	3. 1:n/1 -> strict max-2~3
	4. valid 1:x:2 <-__queue_cut
****/
float __queue_fast(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,
				float *freArr2,float *dbArr2,float *heightArr2,int refLength,
				float light,int *valid,
				int *formatFlag,
				float *fre11,float *fre22,float *fre33,
				float *db11,float *db22,float *db33){

	float fre=0;
	float _fre=0;
	
	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int vs1=0,vs2=0;
	int vk1=0,vk2=0,vk3=0;

	int index1=0;

	if(valid){
		if(*valid&&refLength>5){
			return 0;
		}
	}

	// 1. 3~4 ->1:1 strict max-3~4
	if(length>=3){
		// 012 ->1:1
		for(int i=0;i<length-2;i++){
			if(indexArr[i]+indexArr[i+1]+indexArr[i+2]==3){
				_fre=__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
							&us1,&us2,
							&uk1,&uk2,&uk3);
				
				if(us1==1&&us1==us2){
					if(fabsf(_fre*uk2-freArr[i+1])<5&&
						fabsf(_fre*uk3-freArr[i+2])<5){

						fre=_fre;

						// valid 1:x:2
						if(indexArr[i]==0){
							if(2*uk1==uk3){
								return 0;
							}
						}
					}
				}

				break;
			}
		}

		// 01x ->1:1 
		if(!fre){
			for(int i=0;i<length-2;i++){
				if(indexArr[i]+indexArr[i+1]==1){
					_fre=__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
								&us1,&us2,
								&uk1,&uk2,&uk3);
					
					if(us1==1&&us1==us2){
						if(fabsf(_fre*uk2-freArr[i+1])<5&&
							fabsf(_fre*uk3-freArr[i+2])<5){

							fre=_fre;

							// valid 1:x:2
							if(indexArr[i]==0){
								if(2*uk1==uk3&&refLength>5){
									return 0;
								}
							}
						}
					}

					break;
				}
			}
		}

		// 02x ->1:1 ??? 023
		if(!fre){
			for(int i=0;i<length-2;i++){
				if(indexArr[i]+indexArr[i+1]==2&&
					indexArr[i+2]==3){

					_fre=__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
								&us1,&us2,
								&uk1,&uk2,&uk3);
					
					if(us1==1&&us1==us2){
						if(fabsf(_fre*uk2-freArr[i+1])<5&&
							fabsf(_fre*uk3-freArr[i+2])<5){

							fre=_fre;

							// valid 1:x:2
							if(indexArr[i]==0){
								if(2*uk1==uk3){
									return 0;
								}
							}
						}
					}

					break;
				}
			}
		}

		// 0xx ->1:1
		if(!fre){
			if(indexArr[0]==0){
				_fre=__queue_fre3(freArr[0],freArr[1],freArr[2],
								&us1,&us2,
								&uk1,&uk2,&uk3);
					
				if(us1==1&&us1==us2){
					if(fabsf(_fre*uk2-freArr[1])<5&&
						fabsf(_fre*uk3-freArr[2])<5){

						fre=_fre;

						// valid 1:x:2
						if(2*uk1==uk3){
							return 0;
						}
					}
				}
			}
		}
	}

	// 2. valid-110 ->1-24-5/1-24-7, 124-8-11; 
	// 2.1 valid-80, ->2347 
	if(!fre&&length>=4){
		_fre=__queue_fre3(freArr[0],freArr[1],freArr[2],
				&us1,&us2,
				&uk1,&uk2,&uk3);

		if(uk1){
			__queue_fre3(freArr[1],freArr[2],freArr[3],
						&vs1,&vs2,
						&vk1,&vk2,&vk3);

			if(vk1){
				if(uk1==1&&uk2==2&&uk3==4&&
					((vk1==2&&vk2==4&&vk3==5)||
						(vk1==2&&vk2==4&&vk3==7))){

					fre=_fre;
				}
				else{
					if(uk1==2&&uk2==3&&uk3==4&&
						vk1==3&&vk2==4&&vk3==7){

						fre=freArr[0]/2;
					}
				}
			}
			else{
				if(length>=5){
					if(uk1==1&&uk2==2&&uk3==4){
						__queue_fre3(freArr[2],freArr[3],freArr[4],
									&vs1,&vs2,
									&vk1,&vk2,&vk3);

						if(vk1==4&&vk2==8&&vk3==11){

							fre=freArr[0]/2;
						}
					}
				}
			}

			// valid ->197, 190~202 ❌
			if(!fre&&
				uk1==1&&uk2==2&&uk3==4&&
				dbArr[2]>dbArr[0]&&
				dbArr[2]>dbArr[1]&&
				dbArr[0]>dbArr[1]&&
				freArr[2]>380&&freArr[2]<405){ 

				fre=freArr[2]/2;
			}
		}
	}

	// 3. 1:n/1 -> strict max-2~3
	if(!fre){
		index1=__arr_maxIndex(dbArr, length);

		// valid-6string 80
		if(index1==1&&
			freArr[0]<85&&
			refLength<5){ 

			_fre=__queue_fre2(freArr[0],freArr[1],
							&uk1,&uk2);

			if(uk1==1&&uk2==2&&
				fabsf(freArr[0]*2-freArr[1])<5){

				fre=_fre;
				return fre;
			}
		}

		if(index1+1<length){
			float _fre2=0;

			_fre2=freArr[index1+1];
			_fre=__queue_fre2(freArr[index1],freArr[index1+1],
								&uk1,&uk2);

			if(!uk1&&index1+2<length){ // valid ->110
				_fre2=freArr[index1+2];
				if(_fre2>210&&_fre2<230){
					_fre=__queue_fre2(freArr[index1],freArr[index1+2],
									&uk1,&uk2);
				}
			}

			if(uk1==1&&
				(uk2==2||uk2==3)){
				if(fabsf(_fre*uk2-_fre2)<5){
					fre=_fre;

					if(dbArr[index1]-dbArr[index1+1]>18&&
						freArr[index1]>130){ // 147-17, ->4321string

						return fre;
					}

					if(fre>330&&index1+2<length){ // valid >330
						_fre=__queue_fre3(freArr[index1],freArr[index1+1],freArr[index1+2],
										&us1,&us2,
										&uk1,&uk2,&uk3);
						if(us1){
							if(fabsf(freArr[index1]/uk1*uk2-freArr[index1+1])<5&&
								fabsf(freArr[index1]/uk1*uk3-freArr[index1+2])<5){
								fre=_fre;
							}
						}
					}
				}
			}
		}

		if(!index1){ 
			// valid not 1:2:4/1:3:6
			if(fre&&index1+2<length){
				_fre=__queue_fre3(freArr[index1],freArr[index1+1],freArr[index1+2],
								&us1,&us2,
								&uk1,&uk2,&uk3);

				if(uk1==1&&
					((uk2==2&&uk3==4)||
					(uk2==3&&uk3==6))){ // 124/136 ->__queue_slide, 

					fre=0;
				}
				else if(uk1==1&&uk2==2&&uk3==6&&
					freArr[1]>190&&freArr[1]<204){ // valid ->197 +3,126

					fre=0;
				}
				else{
					if(us1&&us1==2*us2){ // 2:1 ->1:3:4
						if(fabsf(freArr[index1]/uk1*uk2-freArr[index1+1])<5&&
							fabsf(freArr[index1]/uk1*uk3-freArr[index1+2])<5){
							fre=_fre;
						}
					}
				}
			}
		}
		else{
			// valid 1:2:6
			if(fre&&index1==1){
				_fre=__queue_fre3(freArr[0],freArr[1],freArr[2],
								&us1,&us2,
								&uk1,&uk2,&uk3);
				// 126 ->valid
				if(uk1==1&&uk2==2&&uk3==6&&
					_fre<90){

					if(dbArr[1]>dbArr[0]&&
						dbArr[1]-dbArr[2]>20){
						
						if(refLength<6){
							return _fre;
						}
						else{
							if(refLength<8){
								__queue_fre3(freArr2[0],freArr2[1],freArr2[2],
											&us1,&us2,
											&uk1,&uk2,&uk3);

								if(uk1==1&&uk2==2&&uk3==6){
									__queue_fre3(freArr2[1],freArr2[2],freArr2[3],
												&us1,&us2,
												&uk1,&uk2,&uk3);

									if(uk1==2&&uk2==6&&uk3==7){
										return _fre;
									}

									if(freArr2[1]+freArr2[2]<freArr2[3]){
										__queue_fre3(freArr2[1]+freArr2[2],freArr2[3],freArr2[4],
													&us1,&us2,
													&uk1,&uk2,&uk3);

										if(uk1==8&&uk2==11&&uk3==12){
											return _fre;
										}
									}
								}
							}
						}
					}
				}
				// else if(uk1==1&&uk2==2&&uk3==4&&
				// 		_fre<90){ // valid-70~90

				// 	return _fre;
				// }

				__queue_fre3(freArr[1],freArr[2],freArr[3],
							&us1,&us2,
							&uk1,&uk2,&uk3);

				if(uk1==1&&
					((uk2==2&&uk3==4)||
					(uk2==2&&uk3==6))&&
					freArr[2]>190&&freArr[2]<204){ // valid ->197 +3, x126,x124

					fre=0;
				}
			}

			// valid not 1:2:4/1:2:6 1:3:6
			if(fre&&refLength>5){
				for(int i=0;i<index1;i++){
					_fre=__queue_fre2(freArr[i],freArr[index1],
									&uk1,&uk2);

					if(uk1==1&&(uk2==2||uk2==3)){
						_fre=__queue_fre3(freArr[i],freArr[index1],freArr[index1+1],
								&us1,&us2,
								&uk1,&uk2,&uk3);

						if(fabsf(_fre*uk2-freArr[index1])<5&&
							fabsf(_fre*uk3-freArr[index1+1])<5){

							fre=0;
						}

						if(!fre&&index1+2<length){
							_fre=__queue_fre3(freArr[index1],freArr[index1+1],freArr[index1+2],
											&us1,&us2,
											&uk1,&uk2,&uk3);

							if(us1&&us1==3*us2){ // 3:1 ->1:4:5
								if(fabsf(freArr[i]*uk2-freArr[index1+1])<5&&
									fabsf(freArr[i]*uk3-freArr[index1+2])<5){

									fre=_fre;
								}
							}
						}

						// if(_fre<90){ // valid-70~90
						// 	fre=_fre;
						// }

						break;
					}
				}
			}
		}
	}

	// valid ->197
	if(!fre){
		index1=__arr_maxIndex(dbArr, length);

		if(!index1&&
			freArr[index1]>190&&freArr[index1]<204){

			if(dbArr[0]-dbArr[1]>18){
				return freArr[0];
			}
		}
	}

	return fre;
}

/***
	sub<14 ??? 15
	for surper-low fre, and high fre+format nosie, and light-median 
****/
float __queue_direct(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length, float light,int *valid){
	float fre=0;
	float _fre=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int vs1=0,vs2=0;
	int vk1=0,vk2=0,vk3=0;

	int index1=0,index2=0,index3=0,index4=0;

	if(valid){
		if(*valid){
			return 0;
		}
	}

	if(length>=3){
		float *arr1=NULL;
		int _index=0;

		arr1=__vnew(length, NULL);

		_index=__arr_maxIndex(dbArr, length);
		__vsort(dbArr, 3, 1, arr1);

		if(arr1[0]-arr1[2]<12){
			index1=indexArr[0];
			index2=indexArr[1];
			index3=indexArr[2];

			if(index1+index2+index3==3){
				fre=__queue_fre3(freArr[0],freArr[1],freArr[2],
								&us1,&us2,
								&uk1,&uk2,&uk3);
			}

			// test3 [857]; i330 [235] ->1:3:6
			if(fre){ // valid -> filter high fre+format nosie
				int qFlag=0;

				qFlag=__validFre3(freArr[0], freArr[1], freArr[2], fre, uk1, uk2, uk3);
				if(!qFlag){ // 108/207 format noise
					fre=0;

					if(length>=4){
						_fre=__queue_fre3(freArr[1],freArr[2],freArr[3],
										&us1,&us2,
										&uk1,&uk2,&uk3);
						if(uk1==1){
							fre=freArr[1];
						}
					}

					if(!fre){
						if(length>=5){
							_fre=__queue_fre3(freArr[2],freArr[3],freArr[4],
											&us1,&us2,
											&uk1,&uk2,&uk3);
							if(uk1==1){
								fre=freArr[2];
							}
						}
					}
				}
				else{
					if(uk1==1&&uk2==3&&uk3==6){
						if(fre>100){ // high-fre ->193/584/1184
							if(dbArr[0]<dbArr[1]&&
								dbArr[1]<dbArr[2]){

								fre=freArr[1];
							}
							else{
								if(length>=4){
									_fre=__queue_fre3(freArr[1],freArr[2],freArr[3],
													&us1,&us2,
													&uk1,&uk2,&uk3);
									if(uk1==1){
										fre=freArr[1];
									}
								}
							}
						}
					}
					else if(uk1==1&&uk2==2&&uk3==4){ // valid 1:2:4 & 2>4>1
						if(dbArr[1]>dbArr[2]&&dbArr[2]>dbArr[0]){
							// fre=freArr[1];
							return 0;
						}
					}
					else if(2*uk1==uk3){ // valid 1:x:2 ->slide
						return 0;
					}
					else if(2*uk2==uk3){ // x:1:2
						if(length>5){ // ->slide
							fre=0;
						}
						else{
							fre=freArr[1];
						}
					}
					else if(uk1==2&&uk2==4&&uk3==5&&
						_index==2&&
						dbArr[0]<dbArr[1]&&
						freArr[2]>240&&freArr[2]<250){ // ->valid 246

						fre=0;
					}
					else if(uk1==7&&uk2==8&&uk3==12&&
						_index==1&&
						freArr[1]/2>130&&freArr[1]/2<160){ // ->valid 146,x23

						fre=freArr[1]/2;
					}
					else{
						if(uk1!=1&&length>=4){
							_fre=__queue_fre3(freArr[1],freArr[2],freArr[3],
											&us1,&us2,
											&uk1,&uk2,&uk3);
							if(uk1==1){
								fre=freArr[1];
							}
						}
					}
				}
			}
		}

		free(arr1);
	}

	// valid 1:2:4 & 2>4>1
	if(!fre&&length>=3){
		if(indexArr[0]+indexArr[1]+indexArr[2]==3){
			_fre=__queue_fre3(freArr[0],freArr[1],freArr[2],
							&us1,&us2,
							&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==2&&uk3==4){ // valid 1:2:4
				if(dbArr[1]>dbArr[2]&&dbArr[2]>dbArr[0]){
					return 0;
				}
			}
			else if(uk1&&2*uk1==uk3){ // valid 1:x:2 ->slide
				return 0;
			}
		}
	}

	return fre;
}

/***
	for low fre and middle fre, normal ->light medium heavy produce more harmonic
****/
float __queue_slide(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,float light,int *valid,int *status){
	float fre=0;

	int len=0;
	int offset=0;

	int index1=0;
	int k1=0;

	int index2=0;
	int k2=0;

	int jumpFlag=0;

	int tFlag=0;
	int oFlag=0;
	int jFlag=0;

	int c1=0;
	int c2=0;

	if(!length){
		return 0;
	}

	for(int i=0;i<length-2;i++){
		int us1=0,us2=0;
		int uk1=0,uk2=0,uk3=0;

		int _vFlag=1;
		int _index=0;

		index1=0;
		k1=0;

		index2=0;
		k2=0;

		jumpFlag=0;
		len=0;

		_index=__arr_maxIndex(dbArr, length);
		if(_index==i){
			if(dbArr[i]-dbArr[i+1]>18&&
				dbArr[i]-dbArr[i+2]>18){

				_vFlag=0;
			}
		}

		len=__queue_cal(freArr,dbArr,heightArr,length,i,_vFlag,
						&index1,&k1,
						&index2,&k2,
						&jumpFlag);
		if(len){
			__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
						&us1,&us2,
						&uk1,&uk2,&uk3);
			
			if(_index==2&&_index==i+2&&
				freArr[i+2]>220&&(!uk3||uk3==5)&&
				length-2>5){

				i++;
				continue;
			}

			if(len){ // ->329 ,124,max-4
				int _index=0;

				_index=__arr_maxIndex(dbArr, length);
				if(uk1==1&&uk2==2&&uk3==4&&
					_index==i+2&&i>0&&
					freArr[i]>220&&freArr[i]<360){

					return freArr[i+1]/2;
				}	
			}

			if(uk1==2&&uk2==4&&uk3==5&&
				i==0&&
				freArr[2]>240&&freArr[2]<255&&
				length>6){

				__queue_fre3(freArr[2],freArr[3],freArr[4],
						&us1,&us2,
						&uk1,&uk2,&uk3);

				if(uk1==1&&uk2==2&&uk3==3){
					__queue_fre3(freArr[3],freArr[4],freArr[5],
								&us1,&us2,
								&uk1,&uk2,&uk3);

					if(uk1==2&&uk2==3&&uk3==4){

					}

					return freArr[3]/2;
				}
			}
		}

		if(len){
			int f1=0;
			int f2=0;
			float base=0;

			if(len==2){ // twoMove
				tFlag=1;
				fre=__queue_twoMove(freArr,dbArr,heightArr,length,i,
									index1,k1,
									index2,k2,
									jumpFlag,
									&offset);

				if(i==0&&
					fabsf(freArr[0]-fre)<10&&
					jumpFlag==2&&
					dbArr[2]>dbArr[0]&&
					dbArr[2]>dbArr[1]&&
					_index==3){ // 12xn -> xx12

					int k1=0;
					int k2=0;

					__queue_fre2(freArr[2],freArr[3],
								&k1,&k2);

					if(k1==1&&k2==2){
						fre=freArr[3]/2;
					}
				}

				// not k1 & not k2
				if(!fre){ // -> 1:x:2, 1&2 is noise 
					if(length>5&&
						dbArr[i+1]>dbArr[i]&&
						c1<=1){
						
						c1++;
						continue;
					}
				}

				// valid 1:2
				if(!fre){
					if(k1&&k2){
						int _index1=0;
						int _index2=0;

						base=freArr[index2]/k2;
						f1=__queue_has(freArr,length,base,offset,&_index1);
						if(f1){
							fre=base;
						}

						if(!fre){
							base=freArr[index1]/k1;
							f1=__queue_has(freArr,length,base,offset,&_index1);
							if(f1){
								fre=base;
							}
						}

						// if(!fre&&
						// 	i==0&&
						// 	dbArr[1]>dbArr[0]&&
						// 	dbArr[2]>dbArr[0]){

						// 	continue;
						// }

						if(!fre){
							int _index=0;

							if(jumpFlag==1){
								_index=index1+2;
							}
							else{
								_index=index1+1;
							}

							if(dbArr[_index]>dbArr[index1]){
								fre=freArr[_index];
								
								if(i==0&&
									k1==3&&
									_index==2&&
									freArr[2]>238&&freArr[2]<260&&
									dbArr[1]>dbArr[0]&&
									dbArr[2]>dbArr[0]){

									int k1=0;
									int k2=0;

									__queue_fre2(freArr[1], freArr[2], &k1, &k2);
									if(k1==2&&k2==3){
										fre=freArr[1]/2;
									}
								}
							}
						}
					}
				}

				// high-fre 
				if(fre>440){
					fre=freArr[index1]/k1;
				}

				if(fre&&status){
					*status=1;
				}
			}
			else{
				index2=0;
				k2=0;

				if(!jumpFlag){ // oneMove
					oFlag=1;
					fre=__queue_oneMove(freArr,dbArr,heightArr,length,i,
										index1,k1,
										&index2,&k2,
										&offset);

					if(!fre){ // -> x:1:2, no 124/126
						__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
									&us1,&us2,
									&uk1,&uk2,&uk3);
						
						if(length>5&&
							dbArr[i+1]>dbArr[i]&&
							(2*uk2==uk3)&&uk2!=2&&
							c2<=1){ // ||3*uk2==uk3
							
							c2++;
							continue;
						}
						else{
							if(length-i>6&&
								uk1==1&&uk2==2&&(uk3==4||uk3==6)){

								int _index1=0;

								f1=__queue_has(freArr,length,freArr[i],i+1,&_index1);
								if(f1){
									fre=freArr[i];
								}
								else{
									fre=freArr[i+1];
								}
							}
						}
					}
				}
				else{ // jumpMove
					jFlag=1;
					fre=__queue_jumpMove(freArr,dbArr,heightArr,length,i,
										index1,k1,
										jumpFlag,
										&index2,&k2,
										&offset);
				}

				if(!fre){
					float base1=0;
					float base2=0;

					if(k1&&k2){
						int _index1=0;
						int _index2=0;

						int us1=0,us2=0;
						int uk1=0,uk2=0,uk3=0;

						int vs1=0,vs2=0;
						int vk1=0,vk2=0,vk3=0;

						float db1=0;
						float db2=0;

						if(index2-index1>=3){
							int i1=0;
							int i2=0;

							i1=__arr_maxIndex(dbArr+index1, 3)+index1;
							i2=__arr_maxIndex(dbArr+index2, 3)+index2;
							db1=dbArr[i1];
							db2=dbArr[i2];
							if(db1-db2>15){
								fre=freArr[index1]/k1;
							}
						}

						base1=freArr[index1]/k1;
						f1=__queue_has(freArr,length,base1,offset,&_index1);
						if(f1){
							fre=base1;

							if(indexArr[index1]==0&&
								freArr[index1]>120){ // max ->fre=base1

							}
							else{
								f2=__queue_has(freArr,length,freArr[index2]/k2,offset,&_index2);
								if(f2&&_index2<_index1){
									fre=freArr[index2]/k2;
								}
							}
							
							if(fre>440){
								int _k=0;

								_k=util_calRangeTimes(base1, fre, NULL);
								if(_k==2){
									fre=base1;
								}
							}
						}

						if(!fre){ 
							base2=freArr[index2]/k2;
							f1=__queue_has(freArr,length,base2,offset,&_index2);
							if(f1){
								fre=base2;

								if(oFlag){ // -> skip-harmonic
									int _index=0;
									int _k=0;

									_index=__arr_maxIndex(dbArr, length);
									if(index1==0&&_index==0&&
										k1==1&&k2==1){

										_k=util_calRangeTimes(freArr[index1], freArr[index2], NULL);
										if(_k==2||_k==4){
											fre=base1;
										}
									}
								}
							}
						}

						if(!fre){ 
							if(fabsf(base1-base2)<10){	// queue error!!! i31 [107]
								fre=base1;
							}
						}
					}
				}

				if(fre&&status){
					if(oFlag){
						*status=2;
					}
					else{
						*status=3;
					}
				}
			}

			break;
		}
	}

	// if(fre){ // ->skip-harmonic or highFre
	// 	int _index=0;
	// 	int _k=0;

	// 	float _db1=0;
	// 	float _db2=0;

	// 	_index=__arr_maxIndex(dbArr, length);
	// 	if(freArr[_index]<fre){
	// 		_k=util_calRangeTimes(freArr[_index], fre, NULL);
	// 		if(_k==2&&
	// 			fabsf(freArr[_index]*2-fre)<5){
				
	// 			fre=freArr[_index];
	// 		}
	// 	}
	// }

	if(!fre){
		float _fre1=0;
		float _fre2=0;

		float _db1=0;
		float _db2=0;

		int us1=0,us2=0;
		int uk1=0,uk2=0,uk3=0;

		if(k1&&k2){  // priority weak
			_fre1=freArr[index1]/k1;
			_fre2=freArr[index2]/k2;
			
			_db1=dbArr[index1]; // max
			_db2=dbArr[index2]; // min
			
			if(index1==index2){ // 1:2 pk 2:3
				fre=_fre1;


				__queue_fre3(freArr[index1],freArr[index1+1],freArr[index1+2],
							&us1,&us2,
							&uk1,&uk2,&uk3);

				if(uk1==2&&uk2==3&&uk3==4&&
					index1==0&&
					dbArr[0]>dbArr[1]&&
					dbArr[0]>dbArr[2]){

					if(fre>130){ //  260
						fre=_fre1;
					}
					else if(fre>70){ // 30~68? and >=70
						fre=_fre2;
					}
				}
			}
			else{
				if(k1==k2&&index1+1==index2){ // 1:2:4
					_db1=dbArr[index1];
					_db2=dbArr[index2];

					if(_db2-_db1<8){
						fre=_fre1;
					}
					else{
						fre=_fre2;
					}
				}
			}
		}
		else{ // dB
			if(k1){ // only one queue
				int us1=0,us2=0;
				int uk1=0,uk2=0,uk3=0;

				float _fre1=0,_fre2=0,_fre3=0;
				float _db1=0,_db2=0;

				fre=freArr[index1]/k1;

				_fre1=freArr[index1];
				_fre2=freArr[index1+1];
				_fre3=freArr[index1+2];

				_db1=dbArr[index1];
				_db2=dbArr[index1+1];

				if(jumpFlag){
					if(jumpFlag==1){
						_fre2=freArr[index1+2];
						_fre3=freArr[index1+3];

						_db2=dbArr[index1+2];
					}
					else{
						_fre2=freArr[index1+1];
						_fre3=freArr[index1+3];

						_db2=dbArr[index1+1];
					}
				}

				__queue_fre3(_fre1,_fre2,_fre3,
							&us1,&us2,
							&uk1,&uk2,&uk3);

				if(uk1==1&&uk2==2&&(uk3==4||uk3==6)&&
					_db2>_db1){ // 1:2:4/1:2:6

					fre=freArr[index1+1];
				}

				if(fre<50){
					if(dbArr[index1+1]>dbArr[index1]){
						float _fre=0;

						_fre=__queue_fre2(freArr[index1+1],freArr[index1+2],
										&k1,&k2);
						if(k1==1){
							fre=_fre;
						}
					}
				}
			}
		}

		if(fre&&status){
			*status=4;
		}
	}

	// if(!fre&&length>3){
	// 	if(indexArr[0]+indexArr[1]==1){
	// 		__queue_fre2(freArr[0],freArr[1],
	// 					&k1,&k2);

	// 		if(indexArr[0]==0&&k1){
	// 			fre=freArr[0];
	// 		}
	// 	}
	// }

	if(fre&&length>=4){ // ->80, 12-nn,23-nn
		int index1=0;
		int index2=0;

		int k1=0;
		int k2=0;

		index1=__arr_maxIndex(dbArr, length);
		index2=__arr_maxIndex(dbArr+2, length-2)+2;
		if(index1<=1&&
			dbArr[index1]-dbArr[index2]>18&&
			(dbArr[0]-dbArr[index2]>12||
				dbArr[1]-dbArr[index2]>12)){

			float _fre=0;

			_fre=__queue_fre2(freArr[0], freArr[1], &k1, &k2);
			if(((k1==1&&k2==2)||
				(k1==2&&k2==3))&&
				_fre>78&&_fre<85){

				return _fre;
			}
		}
	}

	// valid ->maxDb,  ->skip-harmonic or highFre
	if(fre){
		fre=__queue_slideValid(freArr,dbArr,heightArr,indexArr,length,fre);
	}
	
	// ->197+7 ->190~204
	if(!fre&&length>=8){
		int _index1=0;
		int _index2=0;

		_index1=__arr_maxIndex(dbArr, length);
		if(_index1<3){
			_index2=__arr_maxIndex(dbArr+(_index1+1), length-(_index1+1))+_index1+1;
			if(freArr[_index2]>190&&freArr[_index2]<204){

				return __queue_slide(freArr+_index2, 
									dbArr+_index2, 
									heightArr+_index2, 
									indexArr+_index2, 
									length-_index2, 
									light,valid, status);			
			}
		}
	}

	// valid ->247, valid ->80
	if(fre){
		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if(freArr[_index]>230&&freArr[_index]<255&&
			freArr[_index]/fre>3.6){

			if(_index>0&&
				indexArr[_index-1]==1){

				// __queue_fre2(freArr[_index-1],freArr[_index],
				// 			&k1,&k2);

				// if(k1==2&&k2==3&&
				// 	fabsf(freArr[_index-1]/2*3-freArr[_index])<5){

				// 	if(freArr[_index-1]>140&&freArr[_index-1]<180){ // ->80 

				// 		return freArr[_index-1]/2;
				// 	}
				// }

				if(fabsf(freArr[_index-1]/2-freArr[_index]/3)<5){ // error is 15
					return freArr[_index-1]/2;
				}
			}

			fre=freArr[_index];
		}
	}

	// valid ->110
	if(fre&&
		light>0.98&&
		length>3){

		int _index1=0;
		int _index2=0;

		int k1=0,k2=0;

		_index1=__arr_maxIndex(dbArr, length);
		_index2=__arr_maxIndex(dbArr+2, length-2)+2;
		if(fre>300&&fre<360&&
			_index1==1&&
			_index2==2&&
			dbArr[1]-dbArr[0]<6&&
			dbArr[0]-dbArr[2]>18&&
			fabsf(fre-freArr[_index1])<10){

			__queue_fre2(freArr[0], freArr[1], &k1, &k2);
			if(k1==1&&k2==3){
				fre=freArr[1]/3;
			}
		}
	}

	return fre;
}

static float __queue_slideValid(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,float value){
	float fre=0;
	float fre1=0;

	int flag=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int k1=0;
	int k2=0;

	int _index1=0;
	int _index2=0;

	float *_freArr=NULL;
	float *_dbArr=NULL;
	float *_heightArr=NULL;

	_freArr=__vnew(length, NULL);
	_dbArr=__vnew(length, NULL);
	_heightArr=__vnew(length, NULL);

	memcpy(_freArr, freArr, sizeof(float )*length);
	memcpy(_dbArr, dbArr, sizeof(float )*length);
	memcpy(_heightArr, heightArr, sizeof(float )*length);

	fre=value;
	fre1=fre;

	// dB desc
	__vcorrsort1(_dbArr, _freArr ,_heightArr,NULL,length, 1);

	_index1=__arr_maxIndex(dbArr, length);
	if(fre>freArr[_index1]&&
		fre-freArr[_index1]>10){ // fre not freArr[_index1]
		
		if(_index1==0){
			for(int i=1;i<length;i++){
				if(_freArr[i]>fre-10){
					__queue_fre2(fre,_freArr[i],
								&k1,&k2);

					/***
						110 220 ->660 330
						196 587; 196 (392) 784
					****/
					if(k1){
						if(dbArr[_index1]-_dbArr[i]>10){  // ->110
							flag=1;
							fre=freArr[_index1];
						}
						// else if(k2==2&&
						// 	freArr[_index1]>190&&freArr[_index1]<204&&
						// 	dbArr[_index1]-_dbArr[i]>6){

						// 	flag=__queue_query2(freArr,dbArr,heightArr,length,0,fre/2,fre);
						// 	if(flag){
						// 		fre=fre/2;
						// 	}
						// }
					}

					break;
				}
			}

			if(!flag&&
				_dbArr[0]-_dbArr[1]>24){ // >>

				flag=1;
				fre=freArr[_index1];
			}
		}
		else{
			for(int i=1;i<length;i++){
				if(_freArr[i]>freArr[_index1]-10){
					__queue_fre2(freArr[_index1],_freArr[i],
								&k1,&k2);

					if(k1==1){
						if(dbArr[_index1]-_dbArr[i]>10){
							flag=1;
							fre=freArr[_index1];
						}
						else{ // ->196, 591
							if(_freArr[i]>fre+10&&
								freArr[_index1]>190&&freArr[_index1]<204&&
								dbArr[_index1]-_dbArr[i]>6){

								flag=1;
								fre=freArr[_index1];
							}
						}
					}
					else if(k1==2&&k2==3&&
						fabsf(freArr[_index1]/2*3-_freArr[i])<5&&
						dbArr[_index1]-_dbArr[i]>10){ // 80, 150~170

						flag=1;
						fre=freArr[_index1];
					}

					break;
				}
			}
		}

		if(!flag){
			for(int i=0;i<length-1;i++){
				if(fabsf(fre-freArr[i])<10){ 

					__queue_fre2(freArr[_index1],freArr[i],
								&k1,&k2);

					if(k1==1&&(k2==2||k2==3)){
						if(freArr[_index1]>130){

							if(freArr[_index1]>155&&freArr[_index1]<175&&
								k2==2){ // ->165, filter


							}
							else{
								fre=freArr[_index1];
							}
						}
						else{
							int flag=0;

							if(k2==2){
								flag=__queue_query2(freArr,dbArr,heightArr,length,0,freArr[_index1],fre1,0);
							}
							else{
								flag=__queue_query3(freArr,dbArr,heightArr,length,0,freArr[_index1],fre1,0);
							}

							if(flag){
								fre=freArr[_index1];
							}
						}
					}
					else if(k1==2&&k2==3){ // ->75~90, 190~204
						if((freArr[_index1]>150&&freArr[_index1]<180)||
								(freArr[_index1]>380&&freArr[_index1]<408)){

							fre=freArr[_index1]/2;
						}
					}

					break;
				}
			}
		}

		if(flag&&_index1==0&&
			freArr[0]>100&&freArr[0]<120){

			__queue_fre3(freArr[0],freArr[1],freArr[2],
						&us1,&us2,
						&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==2&&uk3==3){ // ->110

			}
			else{
				flag=__queue_query3(freArr,dbArr,heightArr,length,0,freArr[_index1],fre1,0);
				if(!flag){ // ->330
					fre=value;
				}
			}
		}

		if(!flag&&(_index1==0||_index1==1)&&
			freArr[0]>100&&freArr[0]<120){

			__queue_fre3(freArr[0],freArr[1],freArr[2],
						&us1,&us2,
						&uk1,&uk2,&uk3);

			if(uk1==1&&uk2==2&&uk3==3&&
				fabsf(freArr[0]*2-freArr[1])<5&&
				fabsf(freArr[0]*3-freArr[2])<5){

				flag=1;
				fre=freArr[0];
			}
			else{ // ->{110,220}, valid 100,200; {110,330}; {220,330}
				__queue_fre2(freArr[0],freArr[1],
							&k1,&k2);

				if(k1==1&&k2==2&&
					dbArr[1]-dbArr[2]>18){

					if(length>6){ // >8 valid {100,200}
						flag=__queue_query2(freArr,dbArr,heightArr,length,0,freArr[0],freArr[1],0);
					}
					else{
						flag=1;
					}
					
					if(flag){
						fre=freArr[1]/2;
					}
				}
				else if(length>5){ // valid {220,330}
					__queue_fre2(freArr[1],freArr[2],
								&k1,&k2);

					if(k1==2&&k2==3&&
						dbArr[0]>dbArr[2]&&
						dbArr[1]>dbArr[2]&&
						dbArr[2]-dbArr[3]>20&&
						fabsf(freArr[1]/2*3-freArr[2])<4&&
						fabsf(freArr[0]-freArr[1]/2)<4){

						flag=1;
						fre=freArr[1]/2;
					}
				}
				// else if(k1==1&&k2==3&&
				// 	fabsf(freArr[0]*3-freArr[1])<4){ // valid {110,330}

				// 	flag=__queue_query3(freArr,dbArr,heightArr,length,0,freArr[0],freArr[1]);
				// 	if(flag){
				// 		fre=freArr[1]/3;
				// 	}
				// }
			}
		}
	}

	if(!flag&&
		_index1==0&&
		fre>freArr[_index1]&&
		fre-freArr[_index1]>10&&
		freArr[_index1]>220){ // ->247, 125, refLength<7

		int k1=0,k2=0;

		for(int i=0;i<length;i++){
			if(fabsf(fre-freArr[i])<10&&
				dbArr[0]-dbArr[i]>18){

				__queue_fre2(freArr[0],freArr[i],
							&k1,&k2);

				if(k1==1){
					flag=1;
					fre=freArr[0];
				}

				break;
			}
		}
	}

	// baseDb<2*baseDb||baseDb>2*baseDb
	if(!flag){
		for(int i=0;i<length-1;i++){
			if((indexArr[i]+indexArr[i+1]==1||
					indexArr[i]+indexArr[i+1]==2||
						indexArr[i]+indexArr[i+1]==3)&&
				(fabsf(fre-freArr[i])<10||
					fabsf(fre-freArr[i+1])<10)){ // 01/02/12{03->bug}

				__queue_fre2(freArr[i],freArr[i+1],
							&k1,&k2);

				if(k1==1&&(k2==2||k2==3)){
					if(freArr[i]>130){
						
						if(freArr[i]>155&&freArr[i]<175&&
							k2==2){ // ->165, filter


						}
						else{
							fre=freArr[i];
						}

						break;
					}
					else if(indexArr[i]==0){
						int flag=0;

						if(k2==2){
							flag=__queue_query2(freArr,dbArr,heightArr,length,0,freArr[i],freArr[i+1],0);
						}
						else{
							flag=__queue_query3(freArr,dbArr,heightArr,length,0,freArr[i],freArr[i+1],0);
						}

						if(flag){
							if(fabsf(fre-freArr[i])>10){
								fre=freArr[i];
							}
						}

						break;
					}
				}
				else if(k1==2&&k2==3){
					if(indexArr[i]==0&&
						((freArr[i]>150&&freArr[i]<180)||
							(freArr[i]>380&&freArr[i]<408))){

						fre=freArr[i]/2;
						break;
					}
				}
			}
		}
	}

	free(_freArr);
	free(_dbArr);
	free(_heightArr);

	return fre;
}

/***
	sub<15
	for high fre, light low fre , light for pure
****/
float __queue_weak(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length,float light,int *valid,int *status){
	float fre=0;

	int k1=0;
	int k2=0;

	int k3=0;
	int k4=0;

	if(length<2){
		return 0;
	}

	fre=__queue_weakValid(freArr,dbArr,heightArr,indexArr,length);
	if(fre){
		return fre;
	}

	if(length==2){
		fre=__queue_fre2(freArr[0],freArr[1],
						&k1,&k2);

		if(k1==2&&k2==3){ // 2:3
			if(fabsf(dbArr[0]-dbArr[1])>8){
				fre=0;
			}
		}
		else if(k1==1){
			if(dbArr[0]<dbArr[1]){
				fre=0;
			}
		}
		else{
			if(dbArr[0]>dbArr[1]){
				fre=freArr[0];
			}
			else{
				fre=freArr[1];
			}
		}
	}
	else if(length==3){
		float fre1=0;
		float fre2=0;

		fre1=__queue_fre2(freArr[0],freArr[1],
						&k1,&k2);
		fre2=__queue_fre2(freArr[1],freArr[2],
						&k3,&k4);

		if(k1){
			if(k1==2&&k2==3){ // 2:3
				fre=fre1;

				if(dbArr[0]-dbArr[1]>20&&
					freArr[0]>220){ // 247, 1-max

					fre=freArr[0];
				}

				if(heightArr[0]<5){
					if(dbArr[1]-dbArr[0]>10){
						fre=freArr[1];
					}
					else{
						fre=0;
					}
				}

			}
			else if(k1==1){ // 1:n
				if(k2==2){ // 1:2
					fre=fre2;
					if(fre1<90||dbArr[0]>dbArr[1]){ // valid-70~90
						fre=fre1;
					}
					else if(fre1>300&&
						dbArr[1]-dbArr[0]<2){

						fre=fre1;
					}
				}
				else{
					if(dbArr[1]-dbArr[0]>8){
						fre=fre2;
					}
					else{
						fre=fre1;
					}
				}
			}
		}
		else{ // 1 is noise
			int index=0;

			index=__arr_maxIndex(dbArr, length);
			if(index==0){
				fre=freArr[0];
			}
		}

		if(!fre){
			if(dbArr[0]-dbArr[1]>20&&
				dbArr[1]>dbArr[2]&&
				freArr[0]>220){ // >220

				return freArr[0];
			}
		}

		if(!fre){
			if(k3==1&&k4<4){ // 1:2/1:3
				fre=fre2;
			}
		}

		if(!fre){ // max
			int index=0;

			index=__arr_maxIndex(dbArr, length);
			fre=freArr[index];
		}
	}
	else{
		for(int i=0;i<length-1;i++){
			if(indexArr[i]+indexArr[i+1]==1){ 
				if(fabsf(dbArr[i]-dbArr[i+1])<15){
					fre=__queue_fre2(freArr[i],freArr[i+1],
									&k1,&k2);

					fre=0;
					break;
				}
			}
		}
	}

	if(fre){ // valid 6-times
		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		if(fre>40&&fre<50&&
			freArr[_index]/fre>5.5){

			fre=freArr[_index];
		}
		else if(fre<40&&
			freArr[_index]/fre>7){

			fre=freArr[_index];
		}
	}

	return fre;
}

static float __queue_weakValid(float *freArr,float *dbArr,float *heightArr,int *indexArr,int length){
	float fre=0;

	int index=0;

	int k1=0;
	int k2=0;

	index=__arr_maxIndex(dbArr, length);
	if(index==1){
		if(freArr[index]-freArr[index-1]<50){
			for(int i=index+1;i<length;i++){
				float _fre=0;

				_fre=__queue_fre2(freArr[index],freArr[i],
								&k1,&k2);
				if(k1==1){
					fre=_fre;
					break;
				}
			}
		}
	}

	return fre;
}

/***
	move 12x/1x2
	12x max is 6; 1x2 max is 5
	return 0 fail 1 sucess
****/
static float __queue_twoMove(float *freArr,float *dbArr,float *heightArr,int length,int start,
							int index1,int k1,
							int index2,int k2,
							int jumpFlag,
							int *offset){
	float fre=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int vs1=0,vs2=0;
	int vk1=0,vk2=0,vk3=0;

	// hack ->146, 234; >240, 234/236; ->80, 234
	{
		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		__queue_fre3(freArr[start],freArr[start+1],freArr[start+2],
					&us1,&us2,
					&uk1,&uk2,&uk3);

		if(_index==start&&
			uk1==2&&uk2==3&&uk3==4&&
			freArr[start]/2>130&&freArr[start]/2<160){

			*offset=length-1;
			return freArr[start]/2;
		}
		else if(_index==start&&
			uk1==2&&uk2==3&&(uk3==4||uk3==6)&&
			freArr[start]/2>220&&freArr[start]/2<300){

			*offset=length-1;
			return freArr[start]/2;
		}
		else if(uk1==2&&uk2==3&&
			freArr[start]>150&&freArr[start]<180){

			int flag=0;

			if(dbArr[start+1]>dbArr[start+2]){
				flag=1;
			}
			else{
				flag=__queue_query2(freArr+start, dbArr+start, heightArr+start, length-start, 0, freArr[start]/2, freArr[start],1);
			}

			if(flag){
				return freArr[start]/2;
			}
		}
	}

	*offset=length-1; // reset offset
	for(int i=start+1;i<length-2;i++){
		float _fre1=0;
		float _fre2=0;
		float _fre3=0;

		if(i==start+1||i==start+2||i==start+3){ // 4||5||6
			int f1=0;
			int f2=0;

			__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
						&us1,&us2,
						&uk1,&uk2,&uk3);

			if((i==start+2&&jumpFlag==1)||
				i==start+3){ // compare 1x2/12x

				f1=__queue_isEqual(freArr, length, index1, k1, i, uk1);
				if(f1){
					fre=freArr[index1]/k1; // not freArr[i]/uk1
					break;
				}

				f2=__queue_isEqual(freArr, length, index2, k2, i, uk1);
				if(!f2){
					if(freArr[i]>440||dbArr[i]<dbArr[index2]){
						f2=__queue_isEqual(freArr, length, index2, k2, i, 2*uk1);
					}
				}

				if(f2){
					fre=freArr[index2]/k2; // not freArr[i]/uk1
					break;
				}

				*offset=i+1;
				break;
			}
			else{ // jump
				if(i==start+1){ // 1x2/12x
					_fre1=freArr[i-1];

					if(jumpFlag==1){ // 1x2
						_fre2=freArr[i+1];
					}
					else{ // 12x
						_fre2=freArr[i];
					}

					_fre3=freArr[i+2];
				}
				else{ // start+2 -> 12x
					_fre1=freArr[i-1];
					_fre2=freArr[i+1];
					_fre3=freArr[i+2];
				}
				
				__queue_fre3(_fre1,_fre2,_fre3,
							&vs1,&vs2,
							&vk1,&vk2,&vk3);

				f1=__queue_isEqual(freArr, length, index1, k1, i, uk1);
				f2=__queue_isEqual(freArr, length, index2, k2, i-1, vk1);
				if(!f2){
					if(freArr[i-1]>440||dbArr[i-1]<dbArr[index2]){
						f2=__queue_isEqual(freArr, length, index2, k2, i-1, 2*vk1);
					}
				}

				if(f1||f2){
					int _index1=0;
					int _index2=0;

					if(!(f1&&f2)){ // sucess
						if(f1&&!f2){
							fre=freArr[index1]/k1;
							break;
						}
						else{ // f2&&!f1
							if(vk1==1&&vk2==2&&(vk3==4||vk3==6)){ // 1:2:4/1:2:6
								int _vFlag=0;

								if(i+3<length){
									__queue_fre3(_fre2,_fre3,freArr[i+3],
												&vs1,&vs2,
												&vk1,&vk2,&vk3);

									if(vk1==1&&(vk2==2||vk2==3)){
										_vFlag=__queue_has(freArr,length,_fre2,i+2,&_index1);
									}
								}

								if(_vFlag){ 
									fre=_fre2;

									if(index1==0){ // -> skip-harmonic
										int _index=0;
										int _k=0;

										_index=__arr_maxIndex(dbArr, length);
										if(_index==0&&
											k2==1&&vk1==1){

											_k=util_calRangeTimes(freArr[index1], _fre2, NULL);
											if(_k==2||_k==4){
												fre=freArr[index1];
											}
										}
									}

									break;
								}
							}
							else{
								fre=freArr[index2]/k2;
								if(vs2==4){ // valid ->1:4/1:2:6
									int _vFlag=0;
									float _base=0;

									_base=freArr[index1]/k1;
									_vFlag=__queue_has(freArr,length,_base,i+2,&_index1);
									if(_vFlag){
										fre=_base;
									}
								}

								break;
							}
						}
					}
				}
				else{ // fail
					*offset=i+1;
					break;
				}
			}
		}
	}

	return fre;
}

/***
	one->one
	one->two
	one-jump
****/
static float __queue_oneMove(float *freArr,float *dbArr,float *heightArr,int length,int start,
							int index1,int k1,
							int *index2,int *k2,
							int *offset){
	float fre=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int vs1=0,vs2=0;
	int vk1=0,vk2=0,vk3=0;

	// hack ->240, 236
	{
		int _index=0;

		_index=__arr_maxIndex(dbArr, length);
		__queue_fre3(freArr[start],freArr[start+1],freArr[start+2],
					&us1,&us2,
					&uk1,&uk2,&uk3);

		if(_index==start&&
			uk1==2&&uk2==3&&(uk3==4||uk3==6)&&
			dbArr[start+1]>dbArr[start+2]&&
			freArr[start]/2>220&&freArr[start]/2<300){

			*offset=length-1;
			return freArr[start]/2;
		}
	}

	for(int i=start+1;i<length-2;i++){
		int f1=0;
		int f2=0;

		int _index1=0;
		int _k1=0;

		int _index2=0;
		int _k2=0;

		int _jumpFlag=0;
		int _len=0;

		_len=__queue_cal(freArr,dbArr,heightArr,length, i,0,
						&_index1,&_k1,
						&_index2,&_k2,
						&_jumpFlag);

		*offset=length-1;
		if(_len){
			if(_len==2){ // two
				f1=__queue_isEqual(freArr, length, index1, k1, _index1, _k1);
				if(f1){
					fre=freArr[index1]/k1;
					break;
				}

				*index2=_index1;
				*k2=_k1;

				*offset=i+1;

				break;
			}
			else{
				if(!_jumpFlag){ // one
					f1=__queue_isEqual(freArr, length, index1, k1, _index1, _k1);
					if(f1){ // sucess
						fre=freArr[index1]/k1;
						break;
					}

					*index2=_index1;
					*k2=_k1;

					*offset=i+1;

					break;
				}
				else{ // jump
					fre=__queue_jumpBound(freArr,dbArr,heightArr,length,i,
										_index1,_k1,
										_jumpFlag,
										&_index2,&_k2,
										offset);

					if(!fre){
						f1=__queue_isEqual(freArr, length, index1, k1, _index1, _k1);
						if(f1){ // sucess
							fre=freArr[index1]/k1;
							break;
						}

						if(_k2){
							f2=__queue_isEqual(freArr, length, index1, k1, _index2, _k2);
							if(f2){
								fre=freArr[index1]/k1;
								break;
							}
						}

						*index2=_index1;
						*k2=_k1;

						*offset=_index1+1;

						break;
					}
				}
			}
		}
	}

	return fre;
}

/***
	jump ->one
	jump ->two
	jump ->sucess
****/
static float __queue_jumpMove(float *freArr,float *dbArr,float *heightArr,int length,int start,
							int index1,int k1,
							int jumpFlag,
							int *index2,int *k2,
							int *offset){
	float fre=0;

	fre=__queue_jumpBound(freArr,dbArr,heightArr,length,start,
						index1,k1,
						jumpFlag,
						index2,k2,
						offset);

	if(!fre){
		if(!*k2){ // oneMove
			*index2=0;
			*k2=0;

			fre=__queue_oneMove(freArr,dbArr,heightArr,length,start,
								index1,k1,
								index2,k2,
								offset);
		}
	}

	return fre;
}

/***
	jump ->two
	jump ->one
	jump ->sucess
	12x34: 123+234, 2x3+x34
	1x234: 123+234, x23+234
****/
static float __queue_jumpBound(float *freArr,float *dbArr,float *heightArr,int length,int start,
							int index1,int k1,
							int jumpFlag,
							int *index2,int *k2,
							int *offset){
	float fre=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	int vs1=0,vs2=0;
	int vk1=0,vk2=0,vk3=0;

	*offset=length-1;
	if(start+4<length){ // 5
		float _fre1=0;
		float _fre2=0;
		float _fre3=0;

		int f1=0;
		int f2=0;

		int _index3=0;
		int _uk3=0;

		if(jumpFlag==1){ // 1x2
			_fre1=freArr[start+2];
			_index3=start+2;
		}
		else{ // 12x
			_fre1=freArr[start+1];
			_index3=start+1;
		}

		_fre2=freArr[start+3];
		_fre3=freArr[start+4];

		__queue_fre3(_fre1,_fre2,_fre3,
					&us1,&us2,
					&uk1,&uk2,&uk3);

		_uk3=uk1;

		f1=__queue_isEqual(freArr, length, index1, k1, _index3, uk1);
		
		_fre1=freArr[start+1];
		_fre2=freArr[start+2];
		_fre3=freArr[start+3];
		__queue_fre3(_fre1,_fre2,_fre3,
					&us1,&us2,
					&uk1,&uk2,&uk3);

		_fre1=freArr[start+2];
		_fre2=freArr[start+3];
		_fre3=freArr[start+4];
		__queue_fre3(_fre1,_fre2,_fre3,
					&vs1,&vs2,
					&vk1,&vk2,&vk3);

		f2=__queue_isEqual(freArr, length, start+1, uk1, start+2, uk2);

		if(f1){ 
			if(f2){ //  -> two 
				*index2=start+1;
				*k2=uk1;

				*offset=start+3;
			}
			else{ // -> sucess
				fre=freArr[index1]/k1;

				if(uk1){ // 12x/1x2
					*index2=start+1;
					*k2=uk1;
				}
				else if(vk1&&jumpFlag==2){ // 12x
					*index2=start+2;
					*k2=vk1;
				}

				*offset=start+3;
			}
		}
		else{
			if(f2){ // -> sucess
				fre=freArr[start+1]/uk1;

				*index2=start+1;
				*k2=uk1;

				*offset=start+3;
			}
			else{ 
				if(_uk3||uk1||vk1){ // -> two
					if(_uk3){
						*index2=_index3;
						*k2=_uk3;
					}
					else if(uk1){
						*index2=start+1;
						*k2=uk1;
					}
					else if(vk1){
						*index2=start+2;
						*k2=vk1;
					}

					*offset=*index2+1;
				}
				else{ // -> one
					*offset=start+3;
				}
			}
		}
	}
	else{ // 4 ->index2
		__queue_fre3(freArr[start+1],freArr[start+2],freArr[start+3],
					&us1,&us2,
					&uk1,&uk2,&uk3);
		if(us1){
			*index2=start+1;
			*k2=uk1;
		}
	}

	return fre;
}

static int __queue_isEqual(float *freArr,int length,int index1,int k1,int index2,int k2){
	int flag=0;

	int _index=0;
	int _k=0;

	if(!k1||!k2){
		return 0;
	}

	if(index1==index2){
		if(k1==k2){
			return 1;
		}
		else{
			return 0;
		}
	}

	if(index1>index2){
		_index=index1;
		_k=k1;

		index1=index2;
		k1=k2;

		index2=_index;
		k2=_k;
	}

	_k=util_calRangeTimes(freArr[index1]/k1,freArr[index2],NULL);
	if(_k==k2){
		flag=1;
	}

	return flag;
}

/***
	cal one/two/jump
	jumpFlag 1/2 ->1x2/12x
	return 0 NULL; 1 one; 2 two; 
****/
static int __queue_cal(float *freArr,float *dbArr,float *heightArr,int length,int start,int flag,
					int *index1,int *k1,
					int *index2,int *k2,
					int *jumpFlag){
	int len=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	float base=0;
	int qFlag=0;

	if(start+2>=length){
		return 0;
	}

	base=__queue_fre3(freArr[start],freArr[start+1],freArr[start+2],
				&us1,&us2,
				&uk1,&uk2,&uk3);

	if(flag&&base){
		qFlag=__validFre3(freArr[start],freArr[start+1],freArr[start+2], base,uk1, uk2,uk3);
		if(!qFlag){
			us1=0;
			us2=0;

			uk1=0;
			uk2=0;
			uk3=0;

			if(dbArr[start+2]>dbArr[start+1]&&
				dbArr[start+2]>dbArr[start]){

				return 0;
			}
		}
	}

	if(us1){ // sucess -> one/two 
		len=1;
		if((uk1*2==uk3)||
			(uk1*2==uk2&&uk1!=1)){ // 1:2

			if(uk1*2==uk3){ // 1x2
				*jumpFlag=1;
			}
			else{ // 12x
				*jumpFlag=2;
			}

			*index2=start;
			*k2=1;

			len=2;
		}
		else if(uk1==4&&uk3==6){ // 2x3->456
			*jumpFlag=1;

			*index2=start;
			*k2=2;

			len=2;
		}
		
		*index1=start;
		*k1=uk1;
	}
	else{ // fail -> jump
		if(start+3<length){
			__queue_fre2(freArr[start],freArr[start+1],
						&uk1,&uk2);

			if(uk1&&uk1*2==uk2){ // 12x 
				__queue_fre3(freArr[start],freArr[start+1],freArr[start+3],
						&us1,&us2,
						&uk1,&uk2,&uk3);

				if(us1){
					*index1=start;
					*k1=uk1;

					*jumpFlag=2;
					len=1;
				}
			}
			else{
				__queue_fre2(freArr[start],freArr[start+2],
							&uk1,&uk3);

				if(uk1&&uk1*2==uk3){ // 1x2
					__queue_fre3(freArr[start],freArr[start+2],freArr[start+3],
								&us1,&us2,
								&uk1,&uk2,&uk3);

					if(us1){
						*index1=start;
						*k1=uk1;

						*jumpFlag=1;
						len=1;
					}
				}
			}
		}
	}

	return len;
}

static int __queue_has(float *freArr,int length,float baseFre,int start,int *index){
	int flag=0;

	float fre=0;
	int k1=0;

	int us1=0,us2=0;
	int uk1=0,uk2=0,uk3=0;

	for(int i=start;i<length-2;i++){
		fre=__queue_fre3(freArr[i],freArr[i+1],freArr[i+2],
						&us1,&us2,
						&uk1,&uk2,&uk3);

		if(fre){
			k1=util_calRangeTimes(fre,baseFre,NULL);
			if(k1==1){
				flag=1;
				if(index){
					*index=i;
				}

				break;
			}
		}
	}

	return flag;
}

static void __map_add(QueueMap *map,float *freArr,int length,int index,int k){
	int *indexArr=NULL;
	int *kArr=NULL;

	int *numArr=NULL;

	int indexLength=0;
	int capLength=0;

	int flag=0;
	int k1=0;

	indexArr=map->indexArr;
	kArr=map->kArr;

	numArr=map->numArr;

	indexLength=map->length;
	capLength=map->capLength;

	if(indexLength>=capLength){
		printf("QueueMap is full!!! \n");
		return;
	}

	for(int i=0;i<indexLength;i++){
		int _fre1=0;

		_fre1=freArr[indexArr[i]]/kArr[i];
		k1=util_calRangeTimes(_fre1,freArr[index],NULL);
		if(k1==k){
			numArr[i]++;
			map->length++;

			flag=1;

			break;
		}
	}

	if(!flag){
		indexArr[length]=index;
		kArr[length]=k;

		numArr[length]=1;

		map->length++;
	}
}

static int __map_has(QueueMap *map,float *freArr,int length,int index,int k){
	int _index=-1;

	int *indexArr=NULL;
	int *kArr=NULL;

	int *numArr=NULL;

	int indexLength=0;

	int k1=0;

	indexArr=map->indexArr;
	kArr=map->kArr;

	numArr=map->numArr;

	indexLength=map->length;
	for(int i=0;i<indexLength;i++){
		int _fre1=0;

		_fre1=freArr[indexArr[i]]/kArr[i];
		k1=util_calRangeTimes(_fre1,freArr[index],NULL);
		if(k1==k){
			_index=indexArr[i];
			break;
		}
	}

	return _index;
}

static QueueMap *__createMap(int capLength){
	QueueMap *map=NULL;

	map=(QueueMap *)calloc(1, sizeof(QueueMap ));

	map->indexArr=(int *)calloc(capLength, sizeof(int ));
	map->kArr=(int *)calloc(capLength, sizeof(int ));

	map->numArr=(int *)calloc(capLength, sizeof(int ));

	map->capLength=capLength;

	return map;
}

static void __freeMap(QueueMap *queueMap){

	if(queueMap){
		free(queueMap->indexArr);
		free(queueMap->kArr);

		free(queueMap->numArr);

		free(queueMap);
	}
}

static int __validFre3(float fre1,float fre2,float fre3,float base,int k1,int k2,int k3){
	int flag=1;

	float s1=0;
	float s2=0;

	s1=fabsf(base*k2-fre2);
	s2=fabsf(base*k3-fre3);

	if((s1>5&&fre1<880)||s1>10){
		flag=0;
	}

	if(s2>10){
		flag=0;
	}

	return flag;
}

static int __validFre2(float fre1,float fre2,float base,int k1,int k2){
	int flag=1;

	float s1=0;
	float s2=0;

	if((s1>5&&fre1<880)||s1>10){
		flag=0;
	}

	if((s2>5&&fre2<880)||s2>10){
		flag=0;
	}

	return flag;
}

static int __isValidTimes(int k1,int k2,int k3){
	int flag=0;

	int s1=0;
	int s2=0;

	int _v1=0;

	s1=k2-k1;
	s2=k3-k2;
	if(s1>s2){
		_v1=s1;
		s1=s2;
		s2=_v1;
	}

	if((s1==1&&(s2==1||s2==2||s2==3||s2==4))||
		(s1==2&&(s2==2||s2==3))){

		flag=1;
	}

	return flag;
}

static float __checkFre(float fre1,float fre2,float fre3){
	float fre=0;

	float arr[3]={0,0,0};

	float sub1=0,sub2=0;
	float base=0;

	arr[0]=fre1;
	arr[1]=fre2;
	arr[2]=fre3;

	// asc
	__vsort(arr,3,0,NULL);

	if(arr[2]<90){
		base=0.92;
	}
	else if(arr[2]<120){
		base=1.5;
	}
	else{
		base=2;
	}

	sub1=arr[1]-arr[0];
	sub2=arr[2]-arr[1];

	if(sub1>base||sub2>base||1){
		if(sub1<sub2){
			fre=(arr[0]+arr[1])/2;
		}
		else{
			fre=(arr[1]+arr[2])/2;
		}
	}
	else{
		fre=(arr[0]+arr[1]+arr[2])/3;
	}
	
	return fre;
}

static float __checkFre2(float fre1,float fre2,float fre3){

	return (fre1+fre2+fre3)/3;
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

static int __arr_cut(float *arr,int length,float value){
	int cutLen=0;

	cutLen=length;
	for(int i=0;i<length;i++){
		if(arr[i]>value){
			cutLen=i;
			break;
		}
	}

	return cutLen;
}







