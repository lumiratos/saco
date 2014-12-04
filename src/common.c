/*------------------------------------------------------------------------------

Copyright 2009-2012 Luís M. O. Matos (luismatos@ua.pt), All Rights Reserved.

These programs are supplied free of charge for research purposes only,
and may not be sold or incorporated into any commercial product. There is
ABSOLUTELY NO WARRANTY of any sort, nor any undertaking that they are
fit for ANY PURPOSE WHATSOEVER. Use them at your own risk. If you do
happen to find a bug, or have modifications to suggest, please report
the same to Luis M. O. Matos, luismatos@ua.pt. The copyright notice above
and this statement of conditions must remain an integral part of each
and every copy made of these files.

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "mafImg.h"
#include "mem.h"

/*------------------------------------------------------------------------------
    Pow function from http://martin.ankerl.com/2007/10/04/
      optimized-pow-approximation-for-java-and-c-c/
------------------------------------------------------------------------------*/

double Pow(double a, double b)
	{
	int tmp = (*(1 + (int *)&a));
	int tmp2 = (int)(b * (tmp - 1072632447) + 1072632447);
	double p = 0.0;
	*(1 + (int * )&p) = tmp2;
	return p;
	}
	
/*----------------------------------------------------------------------------*/
	
int BaseToSymbol(int base)
	{
	switch(base)
		{
		case 'a': case 'A': return 0;
		case 'c': case 'C': return 1;
		case 'g': case 'G': return 2;
		case 't': case 'T': return 3;
		default: return 4; //N,n,-, .. //N,n,-, ...
		//case 'n': case 'N': return 4;
		//default: return 5; //N,n,-, .. //N,n,-, ...
		
		}
	}

/*----------------------------------------------------------------------------*/

int BaseTransform(int base)
	{
	switch(base)
		{
		case 'a': case 'A': return 'A';
		case 'c': case 'C': return 'C';
		case 'g': case 'G': return 'G';
		case 't': case 'T': return 'T';
		default: return '-';
		}
	}

/*----------------------------------------------------------------------------*/

int SymbolToBase(int symbol)
	{
	switch(symbol)
		{
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		//case 4: return 'N';
		default: return '-'; //N,n,-, .. //N,n,-, ...
		}
	}

/*----------------------------------------------------------------------------*/

void OutputSymbol(int symbol, FILE *stream)
	{
	switch(symbol)
		{
		case 0: putc('A', stream); break;
		case 1: putc('C', stream); break;
		case 2:	putc('G', stream); break;
		case 3:	putc('T', stream); break;
		default: putc('-', stream); break; // N, n, -, ..
		}
	}
	
/*----------------------------------------------------------------------------*/
	
void StoreSymbol(UChar **mafImgRow, int mafImgRowSize, int base)
	{
	if(mafImgRowSize % ROW_BLOCK_SIZE == 0)
		*mafImgRow = (UChar *)Realloc(*mafImgRow, sizeof(UChar) *
		  (mafImgRowSize + ROW_BLOCK_SIZE), sizeof(UChar) *
		  ROW_BLOCK_SIZE);

	(*mafImgRow)[mafImgRowSize++] = BaseToSymbol(base);
	}
	
/*----------------------------------------------------------------------------*/
void StoreOriginalSymbol(UChar **mafImgRow, int mafImgRowSize, int base)
	{
	if(mafImgRowSize % ROW_BLOCK_SIZE == 0)
		*mafImgRow = (UChar *)Realloc(*mafImgRow, sizeof(UChar) *
	  	(mafImgRowSize + ROW_BLOCK_SIZE), sizeof(UChar) *
	  	ROW_BLOCK_SIZE);

	(*mafImgRow)[mafImgRowSize++] = base;
	}

/*----------------------------------------------------------------------------*/

void AddRowToMAFImg(MAFImg *mafImg, UChar *newRow, int rowSize)
	{
	if(mafImg->nCols && rowSize != mafImg->nCols)
		{
		fprintf(stderr, "Error: trying to add a row with different size\n");
		exit(1);
		}

	mafImg->data = (UChar **)Realloc(mafImg->data, (mafImg->nRows + 1) *
	  sizeof(UChar *), sizeof(UChar *));

	mafImg->data[mafImg->nRows] = newRow;
	mafImg->nRows++;
	if(!mafImg->nCols)
		mafImg->nCols = rowSize;
	}
	
/*----------------------------------------------------------------------------*/

void updateAncestorLine(MAFImg *mafImg, int row, int col)
	{
	int counts[N_SYMBOLS], i;

	// Update only the current column
	// Reset counters
	for(i=0; i < N_SYMBOLS; i++) 
		counts[i] = 0;
		
	// Loop all lines, except the first one which is the ancestorLine
	for(i=1; i <= row; i++) 
		counts[GetMAFPixel(mafImg, i, col)]++;
				
	SetMAFPixel(mafImg, 0, col, 0);
		
	// Set the base which occurs more often
	for(i=1; i < N_SYMBOLS; i++)
		{
		//if(counts[i] > counts[GetMAFPixel(mafImg, 0, col)])
		if(counts[i] >= counts[GetMAFPixel(mafImg, 0, col)])
			SetMAFPixel(mafImg, 0, col, i);
		}
	}
/*----------------------------------------------------------------------------*/
