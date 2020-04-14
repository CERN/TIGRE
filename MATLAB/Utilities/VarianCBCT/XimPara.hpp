#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <iostream>

//typedef struct XimPara

#ifndef STR_XIM
#define STR_XIM
//struct XimPara
typedef struct XimPara
{
	char FileName[256]; 
	int ImgWidth;					// Image Width
	int ImgHeight;					// Image Height
	int PixelNO;

	int BytesPerPixel;				// Determine how to read the data
	int Compression_Indicator;		// Data number in Rec Image Matrix

	double GantryRtn;				// Gantry rotation angle
}XimPara;
#endif

//#ifndef cReadXim_FUN
//#define cReadXim_FUN
// int cReadXim(char *XimFullFile, XimPara *XimStr, int *XimImg);
//#endif
