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
#ifndef MAF_IMG_H_INCLUDED
#define MAF_IMG_H_INCLUDED

typedef unsigned char UChar;

typedef struct
	{
	int		nRows;
	int		nCols;
	UChar	**data;
	}
MAFImg;

typedef struct
	{
	int row;
	int col;
	}
MAFImgCoords;

int GetMAFPixel(MAFImg *mafImg, int row, int col);
void SetMAFPixel(MAFImg *mafImg, int row, int col, int MAFpixel);
void ResetMAFImg(MAFImg *mafImg);
void FreeMAFImg(MAFImg *mafImg);
MAFImg *CreateMAFImg();
MAFImg *CreateMAFImgWith(int rows, int cols);

#endif /* MAF_IMG_H_INCLUDED */