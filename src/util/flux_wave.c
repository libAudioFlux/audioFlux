// 

#include <string.h>
#include <math.h>

#include "../vector/flux_vector.h"
#include "../vector/flux_vectorOp.h"

#include "flux_wave.h"

typedef struct  {
	// RIFF chunk 12
	char riff[4]; // "riff"
	unsigned int size;
	char wav[4];  // "WAVE"

	// fmt sub-chunk >=20
	char fmt[4];  // "fmt "
	unsigned int fmtSize; // >=16 

	unsigned short format; // 1/3/... 
	unsigned short channelNum;

	unsigned int samplate;
	unsigned int byteRate;

	unsigned short blockAlign; 
	unsigned short bit;

	// data sub-chunk 8
	char data[4]; // "data"
	unsigned int dataSize;
} WavHeader;

struct OpaqueWaveRead{
	FILE *fp;

	int samplate;
	int bit;
	int channelNum;

	unsigned short format;

	int dataLength; // head data totalLength
	int curLength; // read
};

struct OpaqueWaveWrite{
	FILE *fp;

	int samplate;
	int bit;
	int channelNum;

	unsigned int fmt_size; 

	int curLength; // write
};

char __wav_header[44]={ 
						0x52, 0x49, 0x46, 0x46, 0x00, 0x00, 0x00, 0x00,
						0x57, 0x41, 0x56, 0x45, 0x66, 0x6d, 0x74, 0x20,
						0x10, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
						0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
						0x00, 0x00, 0x00, 0x00, 0x64, 0x61, 0x74, 0x61,
						0x00, 0x00, 0x00, 0x00 
};

int waveReadObj_new(WaveReadObj *waveObj,char *fileName){
	int status=0;
	WaveReadObj wave=NULL;

	FILE *fp=NULL;
	WavHeader header;

	unsigned int offset=0;
	unsigned int dataSize=0;

	wave=*waveObj=(WaveReadObj )calloc(1, sizeof(struct OpaqueWaveRead ));
	fp=fopen(fileName, "rb+");

	fread(&header, 1, sizeof(header), fp);
	if(header.fmtSize<16){
		printf("wav is error!!!\n");

		fclose(fp);
		free(wave);
		return 1;
	}

	// printf("fmt_size is %d\n",header.fmt_size);
	// printf("format is %d\n",header.format);

	offset=20+header.fmtSize+4; // ???
	fseek(fp, offset, SEEK_SET);
	fread(&dataSize, 1, sizeof(unsigned int ), fp);

	wave->fp=fp;

	wave->samplate=header.samplate;
	wave->channelNum=header.channelNum;
	wave->bit=header.bit;

	wave->format=header.format;

	wave->dataLength=dataSize/(header.bit/ 8);

	return status;
}

int waveReadObj_getInfor(WaveReadObj waveObj,int *samplate,int *bit,int *channelNum){
	int dataLength=0;

	if(samplate){
		*samplate=waveObj->samplate;
	}

	if(channelNum){
		*channelNum=waveObj->channelNum;
	}

	if(bit){
		*bit=waveObj->bit;
	}

	dataLength=waveObj->dataLength;
	return dataLength;
}

int waveReadObj_read(WaveReadObj waveObj,float *dataArr,int dataLength1){
	FILE *fp=NULL;
	int bit=0;

	int channelNum=0;

	unsigned short format;

	int dataLength=0;
	int curLength=0;

	int subLen=0;

	fp=waveObj->fp;
	bit=waveObj->bit;
	if(bit!=8&&bit!=16&&bit!=32){
		return 0;
	}

	channelNum=waveObj->channelNum;
	if(dataLength1%channelNum!=0){
		return 0;
	}

	format=waveObj->format;
	if(format!=1&&format!=3){ // PCM/IEEE float
		return 0;
	}

	dataLength=waveObj->dataLength;
	curLength=waveObj->curLength;

	if(curLength+dataLength1>dataLength){
		subLen=dataLength-(curLength+dataLength1);
	}

	for(int i=curLength,j=0;i<curLength+dataLength1-subLen;i++,j++){
		if(bit==8){
			char v1=0;

			fread(&v1, 1, sizeof(char ), fp);
			dataArr[j]=1.0*v1/(1<<7);
		}
		else if(bit==16){
			short v1=0;

			fread(&v1, 1, sizeof(short ), fp);
			dataArr[j]=1.0*v1/(1<<15);
		}
		else{
			if(format==1){ // 1 --> PCM/uncompressed
				int v1=0;

				// 1<<31 -2147483648
				fread(&v1, 1, sizeof(int ), fp);
				dataArr[j]=-1.0*v1/(1<<31);
			}
			else { // 3 --> IEEE float ???
				fread(dataArr+j, 1, sizeof(float ), fp);
			}
		}
	}

	waveObj->curLength=curLength+dataLength1-subLen;
	return dataLength1-subLen;
}

void waveReadObj_free(WaveReadObj waveObj){

	if(waveObj){
		fclose(waveObj->fp);
		free(waveObj);
	}
}

int waveWriteObj_new(WaveWriteObj *waveObj,char *fileName,
					int *samplate,int *bit,int *channelNum){
	int status=0;
	WaveWriteObj wave=NULL;

	FILE *fp=NULL;
	WavHeader header;

	int _samplate=32000;
	int _bit=16;
	int _channelNum=1;

	wave=*waveObj=(WaveWriteObj )calloc(1, sizeof(struct OpaqueWaveWrite ));
	fp=fopen(fileName, "wb+");

	if(samplate){
		if(*samplate>0&&*samplate<=196000){
			_samplate=*samplate;
		}
	}

	if(bit){
		if(*bit==8||*bit==16||*bit==32){
			_bit=*bit;
		}
	}

	if(_channelNum){
		if(*channelNum>0){
			_channelNum=*channelNum;
		}
	}

	memcpy(&header, __wav_header, sizeof(header));
	header.channelNum=_channelNum;
	header.bit=_bit;
	header.samplate=_samplate;
	header.byteRate = _samplate*_channelNum*(_bit/8);
	header.blockAlign= _channelNum*(_bit/8);

	fwrite(&header, 1, sizeof(header), fp);

	wave->fp=fp;

	wave->samplate=_samplate;
	wave->bit=_bit;
	wave->channelNum=_channelNum;

	wave->samplate=header.samplate;
	wave->channelNum=header.channelNum;
	wave->bit=header.bit;

	wave->fmt_size=header.fmtSize;

	return status;
}

int waveWriteObj_write(WaveWriteObj waveObj,float *dataArr,int dataLength){
	FILE *fp=NULL;

	int bit=0;
	int channelNum=0;

	int curLength=0;

	unsigned int offset=0;
	unsigned int _size=0;
	unsigned int _dataSize=0;

	fp=waveObj->fp;

	bit=waveObj->bit;
	channelNum=waveObj->channelNum;
	if(dataLength%channelNum!=0){
		return 0;
	}

	curLength=waveObj->curLength;

	for(int i=0;i<dataLength;i++){
		if(bit==8){
			char v1=0;

			if(dataArr[i]>=1.0){
				v1=(1<<7)-1;
			}
			else if(dataArr[i]<=-1.0){
				v1=-(1<<7);
			}
			else{
				v1=round(dataArr[i]*(1<<7));
			}
			
			fwrite(&v1, 1, sizeof(char ), fp);
		}
		else if(bit==16){
			short v1=0;

			if(dataArr[i]>=1.0){
				v1=(1<<15)-1;
			}
			else if(dataArr[i]<=-1.0){
				v1=-(1<<15);
			}
			else{
				v1=round(dataArr[i]*(1<<15));
			}

			fwrite(&v1, 1, sizeof(short ), fp);
		}
		else{ // 32 fmt 1
			/***
				int 1<<31 -2147483648 !!!
			****/
			int v1=0;

			if(dataArr[i]>=1.0){
				v1=-((1<<31)+1);
			}
			else if(dataArr[i]<=-1.0){
				v1=(1<<31);
			}
			else {
				v1=-round(dataArr[i]*(1<<31));
			}
			
			fwrite(&v1, 1, sizeof(int ), fp);
			
			// fwrite(dataArr+i, 1, sizeof(float ), fp);
		}
	}

	_dataSize=(curLength+dataLength)*(bit/8);
	_size=sizeof(WavHeader )-8+_dataSize;

	offset=4;
	fseek(fp, offset, SEEK_SET);
	fwrite(&_size, 1, sizeof(unsigned int ), fp);

	offset=20+waveObj->fmt_size+4;
	fseek(fp, offset, SEEK_SET);
	fwrite(&_dataSize, 1, sizeof(unsigned int ), fp);

	fseek(fp, 0, SEEK_END);

	waveObj->curLength=curLength+dataLength; // update
	return dataLength;
}

void waveWriteObj_free(WaveWriteObj waveObj){

	if(waveObj){
		fclose(waveObj->fp);
		free(waveObj);
	}
}










