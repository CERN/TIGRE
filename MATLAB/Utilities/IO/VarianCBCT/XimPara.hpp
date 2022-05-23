#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <iostream>

// Purpose: To fast read .xim files 
// Method: based on ReadXim.m by Fredrik Nordström 2015
// Date: 2017.07
// Author: Yi Du, yi.du@hotmail.com

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
    int KVNormChamber;          // KV norm chamber reading, date: 2022-05-23
}XimPara;
#endif

