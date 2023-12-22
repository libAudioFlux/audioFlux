// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "../util/flux_util.h"

#include "_queue.h"

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










