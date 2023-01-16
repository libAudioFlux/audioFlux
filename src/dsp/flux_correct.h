

#ifndef FLUX_CORRECT_H
#define FLUX_CORRECT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

/***
	频谱校正:频率/幅值/相位 校正都是基于未恢复幅值进行的!!!
	rect/hann/hamm 三种基本窗的fre/amp校正以及amp恢复
	det 谱线修正量 
	value 幅值校正值
****/
void correct_rect(float cur,float left,float right,float *det,float *value);
void correct_hann(float cur,float left,float right,float *det,float *value);
void correct_hamm(float cur,float left,float right,float *det,float *value);

float correct_getRectRecover();
float correct_getHannRecover();
float correct_getHammRecover();

#ifdef __cplusplus
}
#endif

#endif