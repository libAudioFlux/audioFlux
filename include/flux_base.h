

#ifndef FLUX_BASE_H
#define FLUX_BASE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

// window相关
typedef enum{
	Window_Rect=0, 
	Window_Hann,
	Window_Hamm,

	Window_Blackman,
	Window_Kaiser,

	Window_Bartlett,
	Window_Triang,

	Window_Flattop,
	Window_Gauss,

	Window_Blackman_Harris,
	Window_Blackman_Nuttall,
	Window_Bartlett_Hann,

	Window_Bohman,

	Window_Tukey, // tapered cosine

} WindowType;

typedef enum{
	FilterBand_LowPass=0, 
	FilterBand_HighPass,

	FilterBand_BandPass,
	FilterBand_BandStop, // Rejection

	// FilterBand_AllPass,

} FilterBandType;

// spectrum&&spectrogram 相关
typedef enum{
	SpectralData_Power=0,
	SpectralData_Mag,

} SpectralDataType;

typedef enum{
	SpectralFilterBankScale_Linear=0,
	SpectralFilterBankScale_Linspace,

	SpectralFilterBankScale_Mel, 
	SpectralFilterBankScale_Bark,
	SpectralFilterBankScale_Erb, 

	SpectralFilterBankScale_Octave, // similar Constant-Q
	SpectralFilterBankScale_Log,

	SpectralFilterBankScale_Deep, // similar Constant-Q 

	SpectralFilterBankScale_Chroma, // stft-chroma

	SpectralFilterBankScale_LogChroma, // similar cqt-chroma
	SpectralFilterBankScale_DeepChroma, // similar cqt-chroma

} SpectralFilterBankScaleType;

typedef enum{
	SpectralFilterBankStyle_Slaney=0, // Triang
	SpectralFilterBankStyle_ETSI, // Bartlett
	SpectralFilterBankStyle_Gammatone, // gammatone

	SpectralFilterBankStyle_Point,
	SpectralFilterBankStyle_Rect,
	
	SpectralFilterBankStyle_Hann, 
	SpectralFilterBankStyle_Hamm,

	SpectralFilterBankStyle_Blackman, 
	SpectralFilterBankStyle_Bohman,

	SpectralFilterBankStyle_Kaiser,
	SpectralFilterBankStyle_Gauss,

} SpectralFilterBankStyleType;

typedef enum{
	SpectralFilterBankNormal_None=0, // same hight

	SpectralFilterBankNormal_Area, // normal(same hight)/same area
	SpectralFilterBankNormal_BandWidth, 
	
} SpectralFilterBankNormalType;

typedef enum{
	SpectralNoveltyMethod_Sub=0, 

	SpectralNoveltyMethod_Entroy, 
	SpectralNoveltyMethod_KL, 
	SpectralNoveltyMethod_IS, 
	
} SpectralNoveltyMethodType;

typedef enum{
	SpectralNoveltyData_Value=0,
	SpectralNoveltyData_Number,
	
} SpectralNoveltyDataType;

typedef enum{
	ChromaDataNormal_None=0,

	ChromaDataNormal_Max, 
	ChromaDataNormal_Min, 

	ChromaDataNormal_P2, 
	ChromaDataNormal_P1, 
	
} ChromaDataNormalType;

typedef enum{
	CepstralRectify_Log=0,
	CepstralRectify_CubicRoot,

} CepstralRectifyType;

typedef enum{
	CepstralEnergy_Replace=0,
	CepstralEnergy_Append,
	CepstralEnergy_Ignore,

} CepstralEnergyType;

typedef enum{
	PaddingPosition_Center=0,
	PaddingPosition_Right,
	PaddingPosition_Left,

} PaddingPositionType;

typedef enum{
	PaddingMode_Constant=0,
	PaddingMode_Reflect,
	PaddingMode_Wrap, // repeat

} PaddingModeType;

typedef enum{
	WaveletContinue_Morse=0,
	WaveletContinue_Morlet,
	WaveletContinue_Bump,

	WaveletContinue_Paul, 
	WaveletContinue_DOG, // DOG
	WaveletContinue_Mexican, // DOG order=2

	WaveletContinue_Hermit,
	WaveletContinue_Ricker,
	
} WaveletContinueType;

typedef enum{
	WaveletDiscrete_Haar=0,
	WaveletDiscrete_Db, // 2~10/20/30/40
	WaveletDiscrete_Sym, // 2~10/20/30
	WaveletDiscrete_Coif, // 1~5

	WaveletDiscrete_FK, // 4/6/8/14/18/22

	/***
		1.1/1.3/1.5
		2.2/2.4/2.6/2.8
		3.1/3.3/3.5/3.7/3.9
		4.4/5.5/6.8
	****/
	WaveletDiscrete_Bior, 
	WaveletDiscrete_DMey,

} WaveletDiscreteType;



#ifdef __cplusplus
}
#endif

#endif