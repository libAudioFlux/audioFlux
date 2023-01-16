// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"
#include "../vector/flux_complex.h"

#include "./coef/__coef_d.h"
#include "./coef/__coef_r.h"

#include "dwt_filterCoef.h"

typedef enum{
	haar = 0,

	db2,db3,db4,db5,db6,db7,db8,db9,db10,
	db20,db30,db40,

	sym2,sym3,sym4,sym5,sym6,sym7,sym8,sym9,sym10,
	sym20,sym30,

	coif1,coif2,coif3,coif4,coif5,

	fk4,fk6,fk8,fk14,fk18,fk22,

	bior11,bior13,bior15,
	bior22,bior24,bior26,bior28,
	bior31,bior33,bior35,bior37,bior39,
	bior44,bior55,bior68,

	dmey,

	unknow,

} ShortType;

static ShortType __calShortType(WaveletDiscreteType waveletType,int t1,int t2);

/***
	haar =db1
	db 2~10/20/30/40
	sym 2~10/20/30
	coif 1~5
	fk4/6/8/14/18/22
	bior 1.1~1.5/2.2~2.8/3.1~3.9/4.4/5.5/6.8
****/
int dwt_filterCoef(WaveletDiscreteType waveletType,int t1,int t2,int coefType,
				float **loArr,float **hiArr){
	int filterLength=0;
	ShortType st=unknow;

	float *lArr=NULL;
	float *hArr=NULL;

	lArr=__vnew(200, NULL);
	hArr=__vnew(200, NULL);

	st=__calShortType(waveletType, t1, t2);
	switch(st){
		case haar:
			if(!coefType){ // dec
				filterLength=sizeof(__haar_loD)/sizeof(float );
				memcpy(lArr,__haar_loD,sizeof(float )*filterLength);
				memcpy(hArr,__haar_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__haar_loR)/sizeof(float );
				memcpy(lArr,__haar_loR,sizeof(float )*filterLength);
				memcpy(hArr,__haar_hiR,sizeof(float )*filterLength);
			}

			break;
		case db2:
			if(!coefType){ // dec
				filterLength=sizeof(__db2_loD)/sizeof(float );
				memcpy(lArr,__db2_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db2_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db2_loR)/sizeof(float );
				memcpy(lArr,__db2_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db2_hiR,sizeof(float )*filterLength);
			}

			break;
		case db3:
			if(!coefType){ // dec
				filterLength=sizeof(__db3_loD)/sizeof(float );
				memcpy(lArr,__db3_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db3_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db3_loR)/sizeof(float );
				memcpy(lArr,__db3_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db3_hiR,sizeof(float )*filterLength);
			}

			break;
		case db4:
			if(!coefType){ // dec
				filterLength=sizeof(__db4_loD)/sizeof(float );
				memcpy(lArr,__db4_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db4_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db4_loR)/sizeof(float );
				memcpy(lArr,__db4_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db4_hiR,sizeof(float )*filterLength);
			}

			break;
		case db5:
			if(!coefType){ // dec
				filterLength=sizeof(__db5_loD)/sizeof(float );
				memcpy(lArr,__db5_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db5_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db5_loR)/sizeof(float );
				memcpy(lArr,__db5_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db5_hiR,sizeof(float )*filterLength);
			}

			break;
		case db6:
			if(!coefType){ // dec
				filterLength=sizeof(__db6_loD)/sizeof(float );
				memcpy(lArr,__db6_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db6_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db6_loR)/sizeof(float );
				memcpy(lArr,__db6_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db6_hiR,sizeof(float )*filterLength);
			}

			break;
		case db7:
			if(!coefType){ // dec
				filterLength=sizeof(__db7_loD)/sizeof(float );
				memcpy(lArr,__db7_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db7_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db7_loR)/sizeof(float );
				memcpy(lArr,__db7_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db7_hiR,sizeof(float )*filterLength);
			}

			break;
		case db8:
			if(!coefType){ // dec
				filterLength=sizeof(__db8_loD)/sizeof(float );
				memcpy(lArr,__db8_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db8_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db8_loR)/sizeof(float );
				memcpy(lArr,__db8_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db8_hiR,sizeof(float )*filterLength);
			}

			break;
		case db9:
			if(!coefType){ // dec
				filterLength=sizeof(__db9_loD)/sizeof(float );
				memcpy(lArr,__db9_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db9_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db9_loR)/sizeof(float );
				memcpy(lArr,__db9_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db9_hiR,sizeof(float )*filterLength);
			}

			break;
		case db10:
			if(!coefType){ // dec
				filterLength=sizeof(__db10_loD)/sizeof(float );
				memcpy(lArr,__db10_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db10_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db10_loR)/sizeof(float );
				memcpy(lArr,__db10_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db10_hiR,sizeof(float )*filterLength);
			}

			break;
		case db20:
			if(!coefType){ // dec
				filterLength=sizeof(__db20_loD)/sizeof(float );
				memcpy(lArr,__db20_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db20_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db20_loR)/sizeof(float );
				memcpy(lArr,__db20_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db20_hiR,sizeof(float )*filterLength);
			}

			break;
		case db30:
			if(!coefType){ // dec
				filterLength=sizeof(__db30_loD)/sizeof(float );
				memcpy(lArr,__db30_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db30_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db30_loR)/sizeof(float );
				memcpy(lArr,__db30_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db30_hiR,sizeof(float )*filterLength);
			}

			break;
		case db40:
			if(!coefType){ // dec
				filterLength=sizeof(__db40_loD)/sizeof(float );
				memcpy(lArr,__db40_loD,sizeof(float )*filterLength);
				memcpy(hArr,__db40_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__db40_loR)/sizeof(float );
				memcpy(lArr,__db40_loR,sizeof(float )*filterLength);
				memcpy(hArr,__db40_hiR,sizeof(float )*filterLength);
			}

			break;

		case sym2:
			if(!coefType){ // dec
				filterLength=sizeof(__sym2_loD)/sizeof(float );
				memcpy(lArr,__sym2_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym2_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym2_loR)/sizeof(float );
				memcpy(lArr,__sym2_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym2_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym3:
			if(!coefType){ // dec
				filterLength=sizeof(__sym3_loD)/sizeof(float );
				memcpy(lArr,__sym3_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym3_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym3_loR)/sizeof(float );
				memcpy(lArr,__sym3_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym3_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym4:
			if(!coefType){ // dec
				filterLength=sizeof(__sym4_loD)/sizeof(float );
				memcpy(lArr,__sym4_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym4_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym4_loR)/sizeof(float );
				memcpy(lArr,__sym4_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym4_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym5:
			if(!coefType){ // dec
				filterLength=sizeof(__sym5_loD)/sizeof(float );
				memcpy(lArr,__sym5_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym5_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym5_loR)/sizeof(float );
				memcpy(lArr,__sym5_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym5_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym6:
			if(!coefType){ // dec
				filterLength=sizeof(__sym6_loD)/sizeof(float );
				memcpy(lArr,__sym6_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym6_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym6_loR)/sizeof(float );
				memcpy(lArr,__sym6_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym6_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym7:
			if(!coefType){ // dec
				filterLength=sizeof(__sym7_loD)/sizeof(float );
				memcpy(lArr,__sym7_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym7_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym7_loR)/sizeof(float );
				memcpy(lArr,__sym7_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym7_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym8:
			if(!coefType){ // dec
				filterLength=sizeof(__sym8_loD)/sizeof(float );
				memcpy(lArr,__sym8_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym8_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym8_loR)/sizeof(float );
				memcpy(lArr,__sym8_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym8_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym9:
			if(!coefType){ // dec
				filterLength=sizeof(__sym9_loD)/sizeof(float );
				memcpy(lArr,__sym9_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym9_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym9_loR)/sizeof(float );
				memcpy(lArr,__sym9_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym9_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym10:
			if(!coefType){ // dec
				filterLength=sizeof(__sym10_loD)/sizeof(float );
				memcpy(lArr,__sym10_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym10_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym10_loR)/sizeof(float );
				memcpy(lArr,__sym10_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym10_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym20:
			if(!coefType){ // dec
				filterLength=sizeof(__sym20_loD)/sizeof(float );
				memcpy(lArr,__sym20_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym20_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym20_loR)/sizeof(float );
				memcpy(lArr,__sym20_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym20_hiR,sizeof(float )*filterLength);
			}

			break;
		case sym30:
			if(!coefType){ // dec
				filterLength=sizeof(__sym30_loD)/sizeof(float );
				memcpy(lArr,__sym30_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym30_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym30_loR)/sizeof(float );
				memcpy(lArr,__sym30_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym30_hiR,sizeof(float )*filterLength);
			}

			break;

		case coif1:
			if(!coefType){ // dec
				filterLength=sizeof(__coif1_loD)/sizeof(float );
				memcpy(lArr,__coif1_loD,sizeof(float )*filterLength);
				memcpy(hArr,__coif1_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__coif1_loR)/sizeof(float );
				memcpy(lArr,__coif1_loR,sizeof(float )*filterLength);
				memcpy(hArr,__coif1_hiR,sizeof(float )*filterLength);
			}

			break;
		case coif2:
			if(!coefType){ // dec
				filterLength=sizeof(__coif2_loD)/sizeof(float );
				memcpy(lArr,__coif2_loD,sizeof(float )*filterLength);
				memcpy(hArr,__coif2_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__coif2_loR)/sizeof(float );
				memcpy(lArr,__coif2_loR,sizeof(float )*filterLength);
				memcpy(hArr,__coif2_hiR,sizeof(float )*filterLength);
			}

			break;
		case coif3:
			if(!coefType){ // dec
				filterLength=sizeof(__coif3_loD)/sizeof(float );
				memcpy(lArr,__coif3_loD,sizeof(float )*filterLength);
				memcpy(hArr,__coif3_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__coif3_loR)/sizeof(float );
				memcpy(lArr,__coif3_loR,sizeof(float )*filterLength);
				memcpy(hArr,__coif3_hiR,sizeof(float )*filterLength);
			}

			break;
		case coif4:
			if(!coefType){ // dec
				filterLength=sizeof(__coif4_loD)/sizeof(float );
				memcpy(lArr,__coif4_loD,sizeof(float )*filterLength);
				memcpy(hArr,__coif4_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__coif4_loR)/sizeof(float );
				memcpy(lArr,__coif4_loR,sizeof(float )*filterLength);
				memcpy(hArr,__coif4_hiR,sizeof(float )*filterLength);
			}

			break;
		case coif5:
			if(!coefType){ // dec
				filterLength=sizeof(__coif5_loD)/sizeof(float );
				memcpy(lArr,__coif5_loD,sizeof(float )*filterLength);
				memcpy(hArr,__coif5_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__coif5_loR)/sizeof(float );
				memcpy(lArr,__coif5_loR,sizeof(float )*filterLength);
				memcpy(hArr,__coif5_hiR,sizeof(float )*filterLength);
			}

			break;

		case fk4:
			if(!coefType){ // dec
				filterLength=sizeof(__fk4_loD)/sizeof(float );
				memcpy(lArr,__fk4_loD,sizeof(float )*filterLength);
				memcpy(hArr,__fk4_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__fk4_loR)/sizeof(float );
				memcpy(lArr,__fk4_loR,sizeof(float )*filterLength);
				memcpy(hArr,__fk4_hiR,sizeof(float )*filterLength);
			}

			break;
		case fk6:
			if(!coefType){ // dec
				filterLength=sizeof(__fk6_loD)/sizeof(float );
				memcpy(lArr,__fk6_loD,sizeof(float )*filterLength);
				memcpy(hArr,__fk6_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__fk6_loR)/sizeof(float );
				memcpy(lArr,__fk6_loR,sizeof(float )*filterLength);
				memcpy(hArr,__fk6_hiR,sizeof(float )*filterLength);
			}

			break;
		case fk8:
			if(!coefType){ // dec
				filterLength=sizeof(__fk8_loD)/sizeof(float );
				memcpy(lArr,__fk8_loD,sizeof(float )*filterLength);
				memcpy(hArr,__fk8_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__fk8_loR)/sizeof(float );
				memcpy(lArr,__fk8_loR,sizeof(float )*filterLength);
				memcpy(hArr,__fk8_hiR,sizeof(float )*filterLength);
			}

			break;
		case fk14:
			if(!coefType){ // dec
				filterLength=sizeof(__fk14_loD)/sizeof(float );
				memcpy(lArr,__fk14_loD,sizeof(float )*filterLength);
				memcpy(hArr,__fk14_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__fk14_loR)/sizeof(float );
				memcpy(lArr,__fk14_loR,sizeof(float )*filterLength);
				memcpy(hArr,__fk14_hiR,sizeof(float )*filterLength);
			}

			break;
		case fk18:
			if(!coefType){ // dec
				filterLength=sizeof(__fk18_loD)/sizeof(float );
				memcpy(lArr,__fk18_loD,sizeof(float )*filterLength);
				memcpy(hArr,__fk18_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__fk18_loR)/sizeof(float );
				memcpy(lArr,__fk18_loR,sizeof(float )*filterLength);
				memcpy(hArr,__fk18_hiR,sizeof(float )*filterLength);
			}

			break;
		case fk22:
			if(!coefType){ // dec
				filterLength=sizeof(__fk22_loD)/sizeof(float );
				memcpy(lArr,__fk22_loD,sizeof(float )*filterLength);
				memcpy(hArr,__fk22_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__fk22_loR)/sizeof(float );
				memcpy(lArr,__fk22_loR,sizeof(float )*filterLength);
				memcpy(hArr,__fk22_hiR,sizeof(float )*filterLength);
			}

			break;

		case bior11:
			if(!coefType){ // dec
				filterLength=sizeof(__bior11_loD)/sizeof(float );
				memcpy(lArr,__bior11_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior11_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior11_loR)/sizeof(float );
				memcpy(lArr,__bior11_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior11_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior13:
			if(!coefType){ // dec
				filterLength=sizeof(__bior13_loD)/sizeof(float );
				memcpy(lArr,__bior13_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior13_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior13_loR)/sizeof(float );
				memcpy(lArr,__bior13_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior13_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior15:
			if(!coefType){ // dec
				filterLength=sizeof(__bior15_loD)/sizeof(float );
				memcpy(lArr,__bior15_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior15_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior15_loR)/sizeof(float );
				memcpy(lArr,__bior15_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior15_hiR,sizeof(float )*filterLength);
			}

			break;

		case bior22:
			if(!coefType){ // dec
				filterLength=sizeof(__bior22_loD)/sizeof(float );
				memcpy(lArr,__bior22_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior22_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior22_loR)/sizeof(float );
				memcpy(lArr,__bior22_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior22_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior24:
			if(!coefType){ // dec
				filterLength=sizeof(__bior24_loD)/sizeof(float );
				memcpy(lArr,__bior24_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior24_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior24_loR)/sizeof(float );
				memcpy(lArr,__bior24_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior24_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior26:
			if(!coefType){ // dec
				filterLength=sizeof(__bior26_loD)/sizeof(float );
				memcpy(lArr,__bior26_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior26_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior26_loR)/sizeof(float );
				memcpy(lArr,__bior26_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior26_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior28:
			if(!coefType){ // dec
				filterLength=sizeof(__bior28_loD)/sizeof(float );
				memcpy(lArr,__bior28_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior28_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior28_loR)/sizeof(float );
				memcpy(lArr,__bior28_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior28_hiR,sizeof(float )*filterLength);
			}

			break;

		case bior31:
			if(!coefType){ // dec
				filterLength=sizeof(__bior31_loD)/sizeof(float );
				memcpy(lArr,__bior31_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior31_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior31_loR)/sizeof(float );
				memcpy(lArr,__bior31_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior31_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior33:
			if(!coefType){ // dec
				filterLength=sizeof(__bior33_loD)/sizeof(float );
				memcpy(lArr,__bior33_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior33_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior33_loR)/sizeof(float );
				memcpy(lArr,__bior33_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior33_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior35:
			if(!coefType){ // dec
				filterLength=sizeof(__bior35_loD)/sizeof(float );
				memcpy(lArr,__bior35_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior35_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior35_loR)/sizeof(float );
				memcpy(lArr,__bior35_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior35_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior37:
			if(!coefType){ // dec
				filterLength=sizeof(__bior37_loD)/sizeof(float );
				memcpy(lArr,__bior37_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior37_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior37_loR)/sizeof(float );
				memcpy(lArr,__bior37_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior37_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior39:
			if(!coefType){ // dec
				filterLength=sizeof(__bior39_loD)/sizeof(float );
				memcpy(lArr,__bior39_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior39_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior39_loR)/sizeof(float );
				memcpy(lArr,__bior39_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior39_hiR,sizeof(float )*filterLength);
			}

			break;

		case bior44:
			if(!coefType){ // dec
				filterLength=sizeof(__bior44_loD)/sizeof(float );
				memcpy(lArr,__bior44_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior44_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior44_loR)/sizeof(float );
				memcpy(lArr,__bior44_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior44_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior55:
			if(!coefType){ // dec
				filterLength=sizeof(__bior55_loD)/sizeof(float );
				memcpy(lArr,__bior55_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior55_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior55_loR)/sizeof(float );
				memcpy(lArr,__bior55_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior55_hiR,sizeof(float )*filterLength);
			}

			break;
		case bior68:
			if(!coefType){ // dec
				filterLength=sizeof(__bior68_loD)/sizeof(float );
				memcpy(lArr,__bior68_loD,sizeof(float )*filterLength);
				memcpy(hArr,__bior68_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__bior68_loR)/sizeof(float );
				memcpy(lArr,__bior68_loR,sizeof(float )*filterLength);
				memcpy(hArr,__bior68_hiR,sizeof(float )*filterLength);
			}

			break;

		case dmey:
			if(!coefType){ // dec
				filterLength=sizeof(__dmey_loD)/sizeof(float );
				memcpy(lArr,__dmey_loD,sizeof(float )*filterLength);
				memcpy(hArr,__dmey_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__dmey_loR)/sizeof(float );
				memcpy(lArr,__dmey_loR,sizeof(float )*filterLength);
				memcpy(hArr,__dmey_hiR,sizeof(float )*filterLength);
			}

			break;

		default:	// unknown -> sym4
			if(!coefType){ // dec
				filterLength=sizeof(__sym4_loD)/sizeof(float );
				memcpy(lArr,__sym4_loD,sizeof(float )*filterLength);
				memcpy(hArr,__sym4_hiD,sizeof(float )*filterLength);
			}
			else{ // rec
				filterLength=sizeof(__sym4_loR)/sizeof(float );
				memcpy(lArr,__sym4_loR,sizeof(float )*filterLength);
				memcpy(hArr,__sym4_hiR,sizeof(float )*filterLength);
			}
			break;
	}
	
	*loArr=lArr;
	*hiArr=hArr;

	return filterLength;
}

static ShortType __calShortType(WaveletDiscreteType waveletType,int t1,int t2){
	ShortType st=unknow;

	switch(waveletType){
		case WaveletDiscrete_Haar:
			st=haar;
			break;

		case WaveletDiscrete_Db:
			if(t1==2){
				st=db2;
			}
			else if(t1==3){
				st=db3;
			}
			else if(t1==4){
				st=db4;
			}
			else if(t1==5){
				st=db5;
			}
			else if(t1==6){
				st=db6;
			}
			else if(t1==7){
				st=db7;
			}
			else if(t1==8){
				st=db8;
			}
			else if(t1==9){
				st=db9;
			}
			else if(t1==10){
				st=db10;
			}
			else if(t1==20){
				st=db20;
			}
			else if(t1==30){
				st=db30;
			}
			else if(t1==40){
				st=db40;
			}

			break;

		case WaveletDiscrete_Sym:
			if(t1==2){
				st=sym2;
			}
			else if(t1==3){
				st=sym3;
			}
			else if(t1==4){
				st=sym4;
			}
			else if(t1==5){
				st=sym5;
			}
			else if(t1==6){
				st=sym6;
			}
			else if(t1==7){
				st=sym7;
			}
			else if(t1==8){
				st=sym8;
			}
			else if(t1==9){
				st=sym9;
			}
			else if(t1==10){
				st=sym10;
			}
			else if(t1==20){
				st=sym20;
			}
			else if(t1==30){
				st=sym30;
			}

			break;

		case WaveletDiscrete_Coif:
			if(t1==1){
				st=coif1;
			}
			else if(t1==2){
				st=coif2;
			}
			else if(t1==3){
				st=coif3;
			}
			else if(t1==4){
				st=coif4;
			}
			else if(t1==5){
				st=coif5;
			}

			break;

		case WaveletDiscrete_FK:
			if(t1==4){
				st=fk4;
			}
			else if(t1==6){
				st=fk6;
			}
			else if(t1==8){
				st=fk8;
			}
			else if(t1==14){
				st=fk14;
			}
			else if(t1==18){
				st=fk18;
			}
			else if(t1==22){
				st=fk22;
			}

			break;

		case WaveletDiscrete_Bior:
			if(t1==1){
				if(t2==1){
					st=bior11;
				}
				else if(t2==3){
					st=bior13;
				}
				else if(t2==5){
					st=bior15;
				}
			}
			else if(t1==2){
				if(t2==2){
					st=bior22;
				}
				else if(t2==4){
					st=bior24;
				}
				else if(t2==6){
					st=bior26;
				}
				else if(t2==8){
					st=bior28;
				}
			}
			else if(t1==3){
				if(t2==1){
					st=bior31;
				}
				else if(t2==3){
					st=bior33;
				}
				else if(t2==5){
					st=bior35;
				}
				else if(t2==7){
					st=bior37;
				}
				else if(t2==9){
					st=bior39;
				}
			}
			else if(t1==4&&t2==4){
				st=bior44;
			}
			else if(t1==5&&t2==5){
				st=bior55;
			}
			else if(t1==6&&t2==8){
				st=bior68;
			}

			break;

		case WaveletDiscrete_DMey:	
			st=dmey;
			break;
	}

	return st;
}







