/*------------------------------------------------------------------------------

Copyright 2005-2012 Armando J. Pinho (ap@ua.pt), All Rights Reserved.

These programs are supplied free of charge for research purposes only,
and may not be sold or incorporated into any commercial product. There is
ABSOLUTELY NO WARRANTY of any sort, nor any undertaking that they are
fit for ANY PURPOSE WHATSOEVER. Use them at your own risk. If you do
happen to find a bug, or have modifications to suggest, please report
the same to Armando J. Pinho, ap@ua.pt. The copyright notice above
and this statement of conditions must remain an integral part of each
and every copy made of these files.

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "mafImg.h"

/*----------------------------------------------------------------------------*/

int GetMAFPixel(MAFImg *mafImg, int row, int col)
	{
	if(row >= 0 && row < mafImg->nRows && col >= 0 && col < mafImg->nCols)
		return mafImg->data[row][col];
	else
		return 0;
	}

/*----------------------------------------------------------------------------*/

void SetMAFPixel(MAFImg *mafImg, int row, int col, int MAFpixel)
	{
	if(row >= 0 && row < mafImg->nRows && col >= 0 && col < mafImg->nCols)
		mafImg->data[row][col] = MAFpixel;
	}

/*----------------------------------------------------------------------------*/

void ResetMAFImg(MAFImg *mafImg)
	{
	int n;

	for(n = 0 ; n < mafImg->nRows ; n++)
		Free(mafImg->data[n], (((mafImg->nCols - 1) / ROW_BLOCK_SIZE) + 
		  1) * ROW_BLOCK_SIZE * sizeof(UChar));

	Free(mafImg->data, mafImg->nRows * sizeof(UChar *));
	mafImg->data = NULL;
	mafImg->nRows = mafImg->nCols = 0;
	}

/*----------------------------------------------------------------------------*/

void FreeMAFImg(MAFImg *mafImg)
	{
	int n;

	for(n = 0 ; n < mafImg->nRows ; n++)
		Free(mafImg->data[n], mafImg->nCols*sizeof(UChar));

	Free(mafImg->data, mafImg->nRows * sizeof(UChar *));
	Free(mafImg, sizeof(MAFImg));
	}

/*----------------------------------------------------------------------------*/

MAFImg *CreateMAFImg()
	{
	MAFImg *mafImg = (MAFImg *)Calloc(1,sizeof(MAFImg));

	return mafImg;
	}

/*----------------------------------------------------------------------------*/

MAFImg *CreateMAFImgWith(int rows, int cols)
	{
	int i;
	MAFImg *mafImg = (MAFImg *)Calloc(1,sizeof(MAFImg));
	mafImg->data = (UChar **)Calloc(rows, sizeof(UChar *));

	for(i = 0; i < rows; i++)
		mafImg->data[i] = (UChar *)Calloc(cols, sizeof(UChar));

	mafImg->nRows = rows;
	mafImg->nCols = cols;

	return mafImg;
	}

/*----------------------------------------------------------------------------*/
