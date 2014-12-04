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


#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#define DEFAULT_PMODEL_MAX_COUNT	((1 << (sizeof(ACCounter) * 8)) - 1)
#define STORAGE_BITS_FILE_SIZE		32
#define STORAGE_BITS_PMODEL_DELTA_NUM	6
#define STORAGE_BITS_PMODEL_DELTA_DEN	6
#define STORAGE_BITS_PMODEL_MAX_COUNT	30
#define STORAGE_BITS_N_CMODELS		5
#define STORAGE_BITS_TEMPLATE_SIZE	5
#define STORAGE_BITS_TEMPLATE_POSITION	5
#define STORAGE_BITS_TEMPLATE_ID 5
#define STORAGE_BITS_THRESHOLD 16
#define STORAGE_BITS_INVERTED_REPEATS	1
#define STORAGE_BITS_POWER_FUNCTION	1
#define STORAGE_BITS_HSIZE		30
#define STORAGE_BITS_GAMMA		16
#define STORAGE_BITS_LEFT_ANCESTRAL_SIZE 4
#define STORAGE_BITS_RIGHT_ANCESTRAL_SIZE 4
#define STORAGE_BITS_MAF_ROWS	7 //5	// for the multiz28way we have 28 rows
#define STORAGE_BITS_MAF_COLS	24 // 20	
#define STORAGE_BITS_GAMMA 16

#include "mafImg.h"
#include "defs.h"

double Pow(double a, double b);
int BaseToSymbol(int base);
int BaseTransform(int base);
int SymbolToBase(int symbol);
void OutputSymbol(int symbol, FILE *stream);
void StoreSymbol(UChar **mafImgRow, int mafImgRowSize, int base);
void StoreOriginalSymbol(UChar **mafImgRow, int mafImgRowSize, int base);
void AddRowToMAFImg(MAFImg *mafImg, UChar *newRow, int rowSize);
void updateAncestorLine(MAFImg *mafImg, int row, int col);

#endif /* COMMON_H_INCLUDED */
