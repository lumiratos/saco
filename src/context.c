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
#include <string.h>
#include <math.h>
#include "mafImg.h"
#include "context.h"
#include "common.h"
#include "mem.h"
#include "defs.h"

/*----------------------------------------------------------------------------*/

static HCCounters zeroCounters = {0x00, 0x00, 0x00, 0x00};
static HCCounters auxCounters;

/*
static void InitHashTable(CModel *cModel)
	{ 
	cModel->hTable.entries = (Entry **)Calloc(cModel->hTable.size, 
	  sizeof(Entry *));
	cModel->hTable.counters = (HCCounters **)Calloc(cModel->hTable.size, 
	  sizeof(HCCounters *));
	cModel->hTable.entrySize = (unsigned short *)Calloc(cModel->hTable.size, 
	  sizeof(unsigned short));
	cModel->hTable.nUsedEntries = 0;
	cModel->hTable.nUsedKeys = 0;
	}
*/
/*----------------------------------------------------------------------------*/

static void InitArray(CModel *cModel)
	{
	cModel->array.counters = (ACCounter *)Calloc(cModel->nPModels *
	  cModel->nSymbols, sizeof(ACCounter));
	}

/*----------------------------------------------------------------------------*/

static void InitThreshold(CModel *cModel)
        {
        cModel->threshold.idxOcc = (unsigned *)Calloc(cModel->
	  threshold.sizeThreshold, sizeof(unsigned));
        cModel->threshold.symbols = (unsigned char *)Calloc(cModel->
	  threshold.sizeThreshold, sizeof(unsigned char));
        }

/*----------------------------------------------------------------------------*/

static void InsertKey(HashTable *hTable, unsigned hIndex, unsigned key)
	{
	hTable->entries[hIndex] = (Entry *)Realloc(hTable->entries[hIndex],
	  (hTable->entrySize[hIndex] + 1) * sizeof(Entry), sizeof(Entry));

	if(!hTable->entrySize[hIndex])
		hTable->nUsedEntries++;

	hTable->nUsedKeys++;
	hTable->entries[hIndex][hTable->entrySize[hIndex]].key = key;
	hTable->entrySize[hIndex]++;
	}

/*----------------------------------------------------------------------------*/

static void InsertCounters(HashTable *hTable, unsigned hIndex,
  unsigned nHCCounters, unsigned k, unsigned smallCounters)
	{
	hTable->counters[hIndex] = (HCCounters *)Realloc(hTable->
	  counters[hIndex], (nHCCounters + 1) * sizeof(HCCounters), 
	  sizeof(HCCounters));

	if(k < nHCCounters)
		memmove(hTable->counters[hIndex][k + 1], hTable->counters[hIndex][k],
		  (nHCCounters - k) * sizeof(HCCounters));

	hTable->counters[hIndex][k][0] = smallCounters & 0x03;
	hTable->counters[hIndex][k][1] = (smallCounters & (0x03 << 2)) >> 2;
	hTable->counters[hIndex][k][2] = (smallCounters & (0x03 << 4)) >> 4;
	hTable->counters[hIndex][k][3] = (smallCounters & (0x03 << 6)) >> 6;
	}

/*----------------------------------------------------------------------------*/

static HCCounter *GetHCCounters(HashTable *hTable, unsigned key)
	{
	unsigned k = 0, n;
	unsigned hIndex = key % hTable->size; // The hash index

	for(n = 0 ; n < hTable->entrySize[hIndex] ; n++)
		{ 
		if(hTable->entries[hIndex][n].key == key) // If key found 
			{
			switch(hTable->entries[hIndex][n].counters)
				{
				case 0: return hTable->counters[hIndex][k];

				default:
				auxCounters[0] = hTable->entries[hIndex][n].counters & 0x03;
				auxCounters[1] = (hTable->entries[hIndex][n].counters & (0x03 << 2)) >> 2;
				auxCounters[2] = (hTable->entries[hIndex][n].counters & (0x03 << 4)) >> 4;
				auxCounters[3] = (hTable->entries[hIndex][n].counters & (0x03 << 6)) >> 6;
				return auxCounters;
				}
			}

		if(hTable->entries[hIndex][n].counters == 0)
			k++;

		}

	return NULL;
	}

/*----------------------------------------------------------------------------*/

PModel *CreatePModel(unsigned nSymbols)
	{
	PModel *pModel;

	pModel = (PModel *)Malloc(sizeof(PModel));
	pModel->freqs = (unsigned *)Malloc(nSymbols * sizeof(unsigned));

	return pModel;
	}

/*----------------------------------------------------------------------------*/

FloatPModel *CreateFloatPModel(unsigned nSymbols)
	{
	FloatPModel *floatPModel;

	floatPModel = (FloatPModel *)Malloc(sizeof(FloatPModel));
	floatPModel->freqs = (double *)Malloc(nSymbols * sizeof(double));

	return floatPModel;
	}

/*----------------------------------------------------------------------------*/

void UpdateCModelCounter(CModel *cModel, unsigned pModelIdx, unsigned symbol)
	{
	unsigned n;
	ACCounter *aCounters;

	if(cModel->mode == HASH_TABLE_MODE)
		{
		unsigned char smallCounter;
		unsigned i, k = 0;
		unsigned nHCCounters; // The number of HCCounters in this entry
		unsigned hIndex = pModelIdx % cModel->hTable.size; // The hash index

		for(n = 0 ; n < cModel->hTable.entrySize[hIndex] ; n++)
			{
			if(cModel->hTable.entries[hIndex][n].key == pModelIdx) // If key found 
				{
				// If "counters" is zero, then update the "large" counters.
				if(cModel->hTable.entries[hIndex][n].counters == 0)
					{
					if(++cModel->hTable.counters[hIndex][k][symbol] == 255)
						for(i = 0 ; i < cModel->nSymbols ; i++)
							cModel->hTable.counters[hIndex][k][i] >>= 1;

					return;
					}
				
				smallCounter = (cModel->hTable.entries[hIndex][n].counters >>
				  (symbol << 1)) & 0x03;
				/*
				 * If "counters" is non-zero, then this is at least the
				 * second time that this key is generated. Therefore,
				 * if the "small" counter of the symbol if full (i.e.,
				 * is equal to 3), then the "large" counters have to be
				 * inserted into the right position.
				 */
				if(smallCounter == 3)
					{
					nHCCounters = k;
					for(i = n + 1 ; i < cModel->hTable.entrySize[hIndex] ; i++)
						if(cModel->hTable.entries[hIndex][i].counters == 0)
							nHCCounters++;

					InsertCounters(&cModel->hTable, hIndex, nHCCounters, k,
					  cModel->hTable.entries[hIndex][n].counters);
					cModel->hTable.entries[hIndex][n].counters = 0;
					cModel->hTable.counters[hIndex][k][symbol]++;
					return;
					}

				/*
				 * There is still room for incrementing the "small" counter.
				 */
				else
					{
					smallCounter++;
					cModel->hTable.entries[hIndex][n].counters &=
					  ~(0x03 << (symbol << 1));
					cModel->hTable.entries[hIndex][n].counters |=
					  (smallCounter << (symbol << 1));
					return;
					}

				}

			/* Keeps counting the number of HCCounters in this entry */
			if(!cModel->hTable.entries[hIndex][n].counters)
				k++;

			}

		/* If key not found */
		InsertKey(&cModel->hTable, hIndex, pModelIdx);
		cModel->hTable.entries[hIndex][cModel->hTable.entrySize[hIndex] - 1].
		  counters = (0x01 << (symbol << 1));
		}
	else
		{
		aCounters = &cModel->array.counters[pModelIdx * cModel->nSymbols];
		aCounters[symbol]++;
		if(aCounters[symbol] == cModel->maxCount && cModel->maxCount != 0)
			for(n = 0 ; n < cModel->nSymbols ; n++)
				aCounters[n] >>= 1;

		}

	}

/*----------------------------------------------------------------------------*/

void RemoveCModelCounter(CModel *cModel, unsigned pModelIdx, unsigned char symbol)	
	{
	ACCounter *aCounters;

        aCounters = &cModel->array.counters[pModelIdx * cModel->nSymbols];

        if(aCounters[symbol] > 0)
                aCounters[symbol]--;
	else
		{
		printf("Internal error! (RemoveCModelCounter)\n");
		exit(1);
		}
	}

/*----------------------------------------------------------------------------*/

CModel *CreateCModel(unsigned maxCtxSize, unsigned nSymbols, unsigned
  nCtxSymbols, unsigned deltaNum, unsigned deltaDen, unsigned maxCount,
  unsigned hSize, unsigned threshold)

	{
	CModel *cModel;

	if(!(cModel = (CModel *)calloc(1, sizeof(CModel))))
		{
		fprintf(stderr, "Error: in memory allocation\n");
		exit(1);
		}

	if(maxCtxSize > 16)
		{
		fprintf(stderr, "Error: context size cannot be greater than 16\n");
		exit(1);
		}

	cModel->nPModels = (ULL)pow(nCtxSymbols, maxCtxSize);
	cModel->maxCtxSize = maxCtxSize;
	cModel->ctxSize = maxCtxSize;
	cModel->nSymbols = nSymbols;
	cModel->nCtxSymbols = nCtxSymbols;
	cModel->deltaNum = deltaNum;
	cModel->deltaDen = deltaDen;
	cModel->hTable.size = hSize;
	cModel->mode = ARRAY_MODE;
	cModel->maxCount = maxCount;
	cModel->threshold.sizeThreshold = threshold;
	cModel->threshold.indexThreshold = 0;
	
	InitArray(cModel);

	if(threshold != 0)
		InitThreshold(cModel);

	return cModel;
	}

/*----------------------------------------------------------------------------*/

void FreeCModel(CModel *cModel)
{
	Free(cModel->array.counters, cModel->nPModels*cModel->nSymbols*sizeof(ACCounter));
	Free(cModel->threshold.idxOcc, cModel->threshold.sizeThreshold*sizeof(unsigned));
	Free(cModel->threshold.symbols, cModel->threshold.sizeThreshold*sizeof(unsigned char));
	free(cModel);
}

/*----------------------------------------------------------------------------*/

void ComputePModel(CModel *cModel, PModel *pModel, unsigned pModelIdx)

	{
	int symbol;
	ACCounter *aCounters;
	HCCounter *hCounters;

	pModel->sum = 0;
	if(cModel->mode == HASH_TABLE_MODE)
		{
		if(!(hCounters = GetHCCounters(&cModel->hTable, pModelIdx)))
			hCounters = zeroCounters;

		for(symbol = 0 ; symbol < cModel->nSymbols ; symbol++)
			{
			pModel->freqs[symbol] = cModel->deltaNum + cModel->deltaDen *
			  hCounters[symbol];
			pModel->sum += pModel->freqs[symbol];
			}
		}

	else
		{
		aCounters = &cModel->array.counters[pModelIdx * cModel->nSymbols];
		for(symbol = 0 ; symbol < cModel->nSymbols ; symbol++)
			{
			pModel->freqs[symbol] = cModel->deltaNum + cModel->deltaDen *
			  aCounters[symbol];
			pModel->sum += pModel->freqs[symbol];
			}
		}

	}

/*----------------------------------------------------------------------------*/

void ComputeStaticPModel(MAFImg *mafImg, int row, int col, PModel *pModel)
	{
	int n;

	pModel->sum = N_SYMBOLS;
	for(n = 0; n < N_SYMBOLS; n++)
		pModel->freqs[n] = 1;
	
	// If we are in the first line, we only consider the previous base already encoded...
	if(row == 0)	
		{
		pModel->freqs[GetMAFPixel(mafImg, 0, col-1)]++;
		pModel->sum++;
		}
	else // Computing statistics in the current column
		for(n = 0; n < row; n++)
			{
			pModel->freqs[GetMAFPixel(mafImg, n, col)]++;
			pModel->sum++;
			}
	}

/*----------------------------------------------------------------------------*/

void ComputeStaticPModel2(MAFImg *mafImg, int row, int col, PModel *pModel)
	{
	int n;

	pModel->sum = N_SYMBOLS;
	for(n = 0; n < N_SYMBOLS; n++)
		pModel->freqs[n] = 1;

	// If we are in the first line, we only consider the previous base already encoded...
	if(row == 0)	
		{
		pModel->freqs[GetMAFPixel(mafImg, 0, col-1)]++;
		pModel->sum++;
		}
	else // Computing statistics in the current column
		{
		for(n = 0; n < row; n++)
			{
			pModel->freqs[GetMAFPixel(mafImg, n, col-1)]++;
			pModel->freqs[GetMAFPixel(mafImg, n, col)]++;
			pModel->freqs[GetMAFPixel(mafImg, n, col+1)]++;
			pModel->sum += 3;
			}
		pModel->freqs[GetMAFPixel(mafImg, row, col-1)]++;	
		pModel->sum++;
		}
	}

/*----------------------------------------------------------------------------*/
/*
int PModelUsed(CModel *cModel, unsigned pModelIdx)
	{
	unsigned symbol;
	ACCounter *aCounters;

	if(cModel->mode == HASH_TABLE_MODE)
		return (int)GetHCCounters(&(cModel->hTable), pModelIdx);

	else
		{
		aCounters = &(cModel->array.counters[pModelIdx * cModel->nSymbols]);
		for(symbol = 0 ; symbol < cModel->nSymbols ; symbol++)
			if(aCounters[symbol])
				return 1;
		}

	return 0;
	}
*/
/*----------------------------------------------------------------------------*/

double PModelSymbolNats(PModel *pModel, unsigned symbol)
	{
	return log((double)pModel->sum / pModel->freqs[symbol]);
	}

/*----------------------------------------------------------------------------*/

double FractionOfPModelsUsed(CModel *cModel)
	{
	unsigned pModel, symbol, sum, counter = 0;
	ACCounter *aCounters;
	HCCounter *hCounters;

	sum = 0;
	for(pModel = 0 ; pModel < cModel->nPModels ; pModel++)
		{
		if(cModel->mode == HASH_TABLE_MODE)
			{
			hCounters = GetHCCounters(&(cModel->hTable), pModel);
			if(hCounters)
				counter++;

			}
		else
			{
			aCounters = &(cModel->array.counters[pModel * cModel->nSymbols]);
			for(symbol = 0 ; symbol < cModel->nSymbols ; symbol++)
				sum += aCounters[symbol];

			if(sum != 0)
				counter++;
			}
		}

	return (double)counter / cModel->nPModels;
	}

/*----------------------------------------------------------------------------*/

double FractionOfPModelsUsedOnce(CModel *cModel)

	{
	unsigned pModelIdx;
	unsigned symbol, sum, counter = 0;
	ACCounter *aCounters;
	HCCounter *hCounters;

	sum = 0;
	for(pModelIdx = 0 ; pModelIdx < cModel->nPModels ; pModelIdx++)
		{
		if(cModel->mode == HASH_TABLE_MODE)
			{
			hCounters = GetHCCounters(&(cModel->hTable), pModelIdx);

			if(!hCounters)
				continue;

			for(symbol = 0 ; symbol < cModel->nSymbols ; symbol++)
				sum += hCounters[symbol];

			if(sum == 1)
				counter++;

			}
		else
			{
			aCounters = &(cModel->array.counters[pModelIdx * cModel->nSymbols]);
			for(symbol = 0 ; symbol < cModel->nSymbols ; symbol++)
				sum += aCounters[symbol];

			if(sum == 1)
				counter++;
			}
		}

	return (double)counter / cModel->nPModels;
	}

/*----------------------------------------------------------------------------*/

void HashingStats(CModel *cModel)
	{
	unsigned entry, n, k, emptyEntries = 0, nSmallCounters = 0, maxEntrySize = 0;
	ULL possibleKeys;
	double average = (double)cModel->hTable.nUsedKeys / cModel->hTable.size,
	  deviation = 0;

	for(entry = 0 ; entry < cModel->hTable.size ; entry++)
		{
		deviation += fabs(average - cModel->hTable.entrySize[entry]);

		if(!cModel->hTable.entrySize[entry])
			{
			emptyEntries++;
			continue;
			}

		if(cModel->hTable.entrySize[entry] > maxEntrySize)
			maxEntrySize = cModel->hTable.entrySize[entry];

		k = 0;
		for(n = 0 ; n < cModel->hTable.entrySize[entry] ; n++)
			{
			/*
			 * If "counters" is non-zero, then this key did not require
			 * "large" counters
			 */
			if(cModel->hTable.entries[entry][n].counters != 0)
				nSmallCounters++;

			/*
			 * Keeps counting the number of HCCounters in this entry.
			 * counters = 0 means that the corresponding key has been
			 * generated more than once and, therefore, the set of
			 * counters has been allocated.
			 */
			if(cModel->hTable.entries[entry][n].counters == 0)
				k++;

			}

		}

	possibleKeys = powl((double)cModel->nSymbols, (double)cModel->maxCtxSize);

	printf("Hash size ......... %u\n", cModel->hTable.size);
	printf("Used entries ...... %u\n", cModel->hTable.nUsedEntries);
	printf("Ideal entry size .. %.1f\n", average);
	printf("Deviation ......... %.1f\n", deviation / cModel->hTable.size);
	//printf("Used keys ......... %u [%.2f %% of %llu]\n",
	//  cModel->hTable.nUsedKeys, 100.0 * (double)cModel->hTable.nUsedKeys /
	//  possibleKeys, possibleKeys);
	printf("Used keys ......... %u [%.2f %% of %"PRIu64"]\n",
	  cModel->hTable.nUsedKeys, 100.0 * (double)cModel->hTable.nUsedKeys /
	  possibleKeys, (uint64_t)possibleKeys);

	printf("Small counters .... %u\n", nSmallCounters);
	printf("Large counters .... %u\n", cModel->hTable.nUsedKeys-nSmallCounters);
	printf("Max entry size .... %u\n", maxEntrySize);
	}

/*----------------------------------------------------------------------------*/

/*
 *            2
 *          1 X
 */
//static MAFImgCoords template1[] = {{0, -1}, {-1, 0}};
static MAFImgCoords template1[] = {{0, -1}};

/*
 *            3
 *            2
 *          1 X
 */
static MAFImgCoords template2[] = {{0, -1}, {-1, 0}, {-2, 0}};

/*
 *            4
 *            3
 *            2
 *          1 X
 */
//static MAFImgCoords template3[] = {{0, -1}, {-1, 0}, {-2, 0}, {-3, 0}};

/*
 *          3 2 4
 *          1 X
 */
static MAFImgCoords template3[] = {{0, -1}, {-1, 0}, {-1, -1}, {-1, 1}};

/*
 *            2
 *      4 3 1 X
 */
static MAFImgCoords template4[] = {{0, -1}, {-1, 0}, {0, -2}, {0, -3}};

/*
 *            2
 *    5 4 3 1 X
 */
static MAFImgCoords template5[] = {{0, -1}, {-1, 0}, {0, -2}, {0, -3},
	{0, -4}};

/*
 *            
 *    4 3 2 1 X
 */
static MAFImgCoords template6[] = {{0, -1}, {0, -2}, {0, -3}, {0, -4}};

/*	      3
 *            2
 *            1
 *            X
 *
static MAFImgCoords template7[] = {{-1, 0}, {-2, 0}, {-3, 0}};
*/

/*         
 * 9  8  7  6  5
 * 4  3  2  1  X
 */
static MAFImgCoords template7[] = {{0, -1}, {0, -2}, {0, -3}, {0, -4}, {-1, 0}, 
  {-1, -1}, {-1, -2}, {-1, -3}, {-1, -4}};

/*	        7 3
 *          6 2
 *          5 1
 *          4 X
 */
//static MAFImgCoords template8[] = {{-1, 0}, {-2, 0}, {-3, 0}, 
//	{0, -1}, {-1, -1}, {-2, -1}, {-3,-1}};
	
/*	          4
 *          8 3
 *          7 2
 *          6 1
 *          5 X
 */
static MAFImgCoords template8[] = {{-1, 0}, {-2, 0}, {-3, 0}, {-4, 0}, {0, -1}, 
  {-1, -1}, {-2, -1}, {-3, -1}};


/*	     11 7 3
 *       10 6 2
 *        9 5 1
 *        8 4 X
 */
//static MAFImgCoords template9[] = {{-1, 0}, {-2, 0}, {-3, 0}, 
//	{0, -1}, {-1, -1}, {-2, -1}, {-3, -1}, {0, -2}, 
//	{-1, -2}, {-2, -2}, {-3,-2}};
static MAFImgCoords template9[] = {{-1, 0}, {-2, 0}, {-3, 0}, 
	{0, -1}, {-1, -1}, {-2, -1}, {-3, -1}, {0, -2}, 
	{-1, -2}};

/*
 * From cod:martins:98a, 2-norm (reduced).
 *
 *              6
 *        8  4  2  3  7
 *        5  1  X
 */
static MAFImgCoords template10[] = {{0, -1}, {-1, 0}, {-1, 1}, {-1, -1},
	{0, -2}, {-2, 0}, {-1, 2}, {-1, -2}};

/*
 * From cod:martins:98a, 2-norm (reduced).
 *
 *          10  6  9   
 *        8  4  2  3  7
 *        5  1  X
 */
static MAFImgCoords template11[] = {{0, -1}, {-1, 0}, {-1, 1}, {-1, -1},
	{0, -2}, {-2, 0}, {-1, 2}, {-1, -2}, {-2, 1}, {-2, -1}};

/*
 * From cod:martins:98a, 2-norm (reduced).
 *
 *       12 10  6  9 11
 *        8  4  2  3  7
 *        5  1  X
 */
static MAFImgCoords template12[] = {{0, -1}, {-1, 0}, {-1, 1}, {-1, -1},
	{0, -2}, {-2, 0}, {-1, 2}, {-1, -2}, {-2, 1}, {-2, -1}, {-2, 2},
	{-2, -2}};

/*
 * From cod:martins:98a, 2-norm.
 *
 *             14
 *       12 10  6  9 11
 *        8  4  2  3  7
 *    13  5  1  X
 */
static MAFImgCoords template13[] = {{0, -1}, {-1, 0}, {-1, 1}, {-1, -1},
	{0, -2}, {-2, 0}, {-1, 2}, {-1, -2}, {-2, 1}, {-2, -1}, {-2, 2},
	{-2, -2}, {0, -3}, {-3, 0}};

/*
 * From cod:martins:98a, 2-norm.
 *
 *             14
 *       12 10  6  9 11
 *    16  8  4  2  3  7 15
 *    13  5  1  X
 */
static MAFImgCoords template14[] = {{0, -1}, {-1, 0}, {-1, 1}, {-1, -1},
	{0, -2}, {-2, 0}, {-1, 2}, {-1, -2}, {-2, 1}, {-2, -1}, {-2, 2},
	{-2, -2}, {0, -3}, {-3, 0}, {-1, 3}, {-1, -3}};

/*
 * Non-causal template
 *
 *        2
 *     1  X  3
 *        4
 */
static MAFImgCoords template20[] = {{0, -1}, {-1, 0}, {0, 1}, {1, 0}};

/*
 * Non-causal template
 *
 *       5  2  6
 *       1  X  3
 *       8  4  7
 */
static MAFImgCoords template21[] = {{0, -1}, {-1, 0}, {0, 1}, {1, 0},
    {-1, -1}, {-1, 1}, {1, 1}, {1, -1}};

/*
 * Non-causal template
 *
 *         10
 *       5  2  6
 *    9  1  X  3 11
 *       8  4  7
 *         12
 */
static MAFImgCoords template22[] = {{0, -1}, {-1, 0}, {0, 1}, {1, 0},
    {-1, -1}, {-1, 1}, {1, 1}, {1, -1}, {0, -2}, {-2, 0}, {0, 2}, {2, 0}};

/*
 *           1
 *           X
 */
static MAFImgCoords template23[] = {{-1, 0}};

/*
 *           
 *         1 X
 */
static MAFImgCoords template24[] = {{0, -1}};


/*
 * From cod:martins:98a, 1-norm.
 *
 *             12
 *          11  6 10 16
 *    15  9  5  2  4  8 14
 * 13  7  3  1  X
 */
//static MAFImgCoords template11[] = {{0, -1}, {-1, 0}, {0, -2}, {-1, 1},
//	{-1, -1}, {-2, 0}, {0, -3}, {-1, 2}, {-1, -2}, {-2, 1}, {-2, -1},
//	{-3, 0}, {0, -4}, {-1, 3}, {-1, -3}, {-2, 2}};

/*----------------------------------------------------------------------------*/

int GetPModelIdx(MAFImg *mafImg, int row, int col, CModel *cModel,
  CTemplate *cTemplate)
	{
	int n, mafPixel, idx = 0, prod = 1;

	// Note that the real size used for context calculation may be smaller
	// than the template size. By changing the value of "cModel->ctxSize"
	// this can be dynamically controlled.
	for(n = 0 ; n < cModel->ctxSize ; n++)
		{
		mafPixel = GetMAFPixel(mafImg, row + cTemplate->position[n].row, 
		  col + cTemplate->position[n].col);
		idx += mafPixel * prod;
		prod *= cModel->nSymbols;
		}

	return idx;
	}

/*----------------------------------------------------------------------------*/

int GetPModelIdx2(MAFImg *mafImg, int row, int col, CModel *cModel,
  CTemplate *cTemplate, char alm)
	{
	int n, mafPixel, idx = 0, prod = 1;

	// Note that the real size used for context calculation may be smaller
	// than the template size. By changing the value of "cModel->ctxSize"
	// this can be dynamically controlled.
	for(n = 0 ; n < cModel->ctxSize ; n++)
		{	
		mafPixel = GetMAFPixel(mafImg, row + cTemplate->position[n].row, 
		  col + cTemplate->position[n].col);

		if(alm == 'n' && row + cTemplate->position[n].row == 0)
			mafPixel = 0;

		idx += mafPixel * prod;
		prod *= cModel->nSymbols;
		}

	return idx;
	}

/*----------------------------------------------------------------------------*/

int GetPModelIdx3(MAFImg *mafImg, int row, int col)
	{
	int i, j, idx = 0, prod = 1, n;
	unsigned *count = NULL, *ind = NULL, max = 0;
	
	count = (unsigned *)Calloc(N_SYMBOLS, sizeof(unsigned));
	ind = (unsigned *)Calloc(N_SYMBOLS, sizeof(unsigned));
	
	// Count the symbols in the current column
	for(i=0; i < row; i++)
		count[GetMAFPixel(mafImg, i, col)]++;

	for(i = 0; i < N_SYMBOLS; i++)
		ind[i] = i;
	
	for(i = 0; i < N_SYMBOLS; i++)
		{
		max = count[ind[i]];
		for(j = i+1; j < N_SYMBOLS; j++)
			{
			if(count[ind[j]] > max)
				{
				max = count[ind[j]];
				n = ind[i];
				ind[i] = ind[j];
				ind[j] = n;
				}
			}
			idx += ind[i] * prod;
			prod *= N_SYMBOLS;
		}
	
	free(count);
	free(ind);

	return idx;
	}

/*----------------------------------------------------------------------------*/

CTemplate *InitTemplate(int templateId)
	{
	CTemplate *cTemplate;

	cTemplate = Malloc(sizeof(CTemplate));

	switch(templateId)
		{
		case 1:
			cTemplate->position = template1;
			cTemplate->size = sizeof(template1) / sizeof(template1[0]);
			break;

		case 2:
			cTemplate->position = template2;
			cTemplate->size = sizeof(template2) / sizeof(template2[0]);
			break;

		case 3:
			cTemplate->position = template3;
			cTemplate->size = sizeof(template3) / sizeof(template3[0]);
			break;

		case 4:
			cTemplate->position = template4;
			cTemplate->size = sizeof(template4) / sizeof(template4[0]);
			break;

		case 5:
			cTemplate->position = template5;
			cTemplate->size = sizeof(template5) / sizeof(template5[0]);
			break;

		case 6:
			cTemplate->position = template6;
			cTemplate->size = sizeof(template6) / sizeof(template6[0]);
			break;

		case 7:
			cTemplate->position = template7;
			cTemplate->size = sizeof(template7) / sizeof(template7[0]);
			break;


		case 8:
			cTemplate->position = template8;
			cTemplate->size = sizeof(template8) / sizeof(template8[0]);
			break;

		case 9:
			cTemplate->position = template9;
			cTemplate->size = sizeof(template9) / sizeof(template9[0]);
			break;

		case 10:
			cTemplate->position = template10;
			cTemplate->size = sizeof(template10) / sizeof(template10[0]);
			break;

		case 11:
			cTemplate->position = template11;
			cTemplate->size = sizeof(template11) / sizeof(template11[0]);
			break;

		case 12:
			cTemplate->position = template12;
			cTemplate->size = sizeof(template12) / sizeof(template12[0]);
			break;

		case 13:
			cTemplate->position = template13;
			cTemplate->size = sizeof(template13) / sizeof(template13[0]);
			break;

		case 14:
			cTemplate->position = template14;
			cTemplate->size = sizeof(template14) / sizeof(template14[0]);
			break;

		case 20:
			cTemplate->position = template20;
			cTemplate->size = sizeof(template20) / sizeof(template20[0]);
			break;

		case 21:
			cTemplate->position = template21;
			cTemplate->size = sizeof(template21) / sizeof(template21[0]);
			break;

		case 22:
			cTemplate->position = template22;
			cTemplate->size = sizeof(template22) / sizeof(template22[0]);
			break;

		case 23: 
			cTemplate->position = template23;
                        cTemplate->size = sizeof(template23) / sizeof(template23[0]);
                        break;

		case 24:
                        cTemplate->position = template24;
                        cTemplate->size = sizeof(template24) / sizeof(template24[0]);
                        break;

		default:
			fprintf(stderr, "Error: invalid template id\n");
			exit(1);
		}

	return cTemplate;
	}

/*----------------------------------------------------------------------------*/

CTemplate *InitAncestorTemplate(int leftSize, int rightSize)
	{
	CTemplate *cTemplate;
	int i, c;

	cTemplate = Malloc(sizeof(CTemplate));

	// Allocate memory for the ancestor template line
	cTemplate->position = (MAFImgCoords *)Calloc(abs(leftSize) + 
	  abs(rightSize)+1, sizeof(MAFImgCoords));

	c = 0;
	for(i=1; i <= (leftSize > rightSize ? leftSize : rightSize); i++)
		{
		if(i <= leftSize)
			cTemplate->position[c++].col = -i;

		if(i <= rightSize)
			cTemplate->position[c++].col = +i;
		}
	
	cTemplate->size = abs(leftSize) + abs(rightSize) + 1;

	return cTemplate;
	}

/*----------------------------------------------------------------------------*/

CTemplate *InitAncestorTemplate2(int leftSize, int rightSize)
	{
	CTemplate *cTemplate;
	int i, c;
	
	cTemplate = Malloc(sizeof(CTemplate));
	// Allocate memory for the ancestor template line 
	cTemplate->position = (MAFImgCoords *)Calloc(abs(leftSize) + 
	  abs(rightSize)+1, sizeof(MAFImgCoords));
	
	c = -abs(leftSize);
	for(i = 0; i < abs(leftSize)+abs(rightSize)+1; i++)
		cTemplate->position[i].col = c++;
	
	cTemplate->size = abs(leftSize) + abs(rightSize) + 1;

	return cTemplate;
	}

/*----------------------------------------------------------------------------*/

void ShowTemplate(CTemplate *cTemplate)
	{
	int minRow, maxRow, minCol, maxCol, n, row, col;
	int **templateMatrix;

	minRow = maxRow = cTemplate->position[0].row;
	minCol = maxCol = cTemplate->position[0].col;

	for(n = 1 ; n < cTemplate->size ; n++)
		{
		if(cTemplate->position[n].row > maxRow)
			maxRow = cTemplate->position[n].row;

		if(cTemplate->position[n].row < minRow)
			minRow = cTemplate->position[n].row;

		if(cTemplate->position[n].col > maxCol)
			maxCol = cTemplate->position[n].col;

		if(cTemplate->position[n].col < minCol)
			minCol = cTemplate->position[n].col;
		}

	templateMatrix = (int **)Calloc(maxRow - minRow + 2, sizeof(int *));

	for(row = 0 ; row < maxRow - minRow + 2 ; row++)
		templateMatrix[row] = (int *)Calloc(maxCol - minCol + 2, sizeof(int));
	
	for(n = 0 ; n < cTemplate->size ; n++)
		templateMatrix[cTemplate->position[n].row - minRow]
		  [cTemplate->position[n].col - minCol] = n + 1;

	templateMatrix[-minRow][-minCol] = -1;
	for(row = 0 ; row < maxRow - minRow + 2 ; row++)
		{
		for(col = 0 ; col < maxCol - minCol + 2 ; col++)
			if(templateMatrix[row][col])
				{
				if(templateMatrix[row][col] == -1)
					printf("  X");
				else
					printf("%3d", templateMatrix[row][col]);
				}
			else
				printf("   ");

		putchar('\n');
		}
	
	for(row = 0 ; row < maxRow - minRow + 2 ; row++)
		free(templateMatrix[row]);

	free(templateMatrix);
	}
