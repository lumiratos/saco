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

/*
 * Data structures for handling finite-context models.
 */

#ifndef CONTEXT_H_INCLUDED
#define CONTEXT_H_INCLUDED

#define MAX_ARRAY_MEMORY			1024 /* DNA: size 12 = 128 MB */
#define ARRAY_MODE				0
#define HASH_TABLE_MODE				1

#pragma pack(1)

typedef unsigned short ACCounter; /* Size of context counters for arrays */
typedef unsigned char  HCCounter; /* Size of context counters for hash tables */
typedef unsigned long long ULL;
typedef HCCounter HCCounters[4];

typedef struct
	{
	unsigned	key;		/* The key stored in this entry */
	unsigned char	counters;	/* "Small" counters: 2 bits for each one */
	}
Entry;

typedef struct
	{
	unsigned	size;		/* Size of the hash table */
	unsigned short	*entrySize;	/* Number of keys in this entry */
	Entry		**entries;	/* The heads of the hash table lists */
	HCCounters	**counters;	/* The context counters */
	unsigned	nUsedEntries;
	unsigned	nUsedKeys;
	}
HashTable;

typedef struct
	{
	ACCounter	*counters;
	}
Array;

typedef struct
	{
	unsigned 	sizeThreshold;	/* Threshold size */
	unsigned	indexThreshold; /* Current position in the threshold */
	unsigned 	*idxOcc;	/* Index of the string by number of occurrences */
	unsigned char	*symbols;	/* symbol for each threshold occureence */
	}
Threshold;

typedef struct
	{
	unsigned	*freqs;
	unsigned	sum;
	}
PModel;

typedef struct
	{
	double		*freqs;
	}
FloatPModel;

typedef struct
	{
	unsigned	maxCtxSize;	/* Maximum depth of context template */
	unsigned	ctxSize;	/* Current depth of context template */
	unsigned	nSymbols;	/* Number of coding symbols */
	unsigned	nCtxSymbols;	/* Number of symbols used for context computation */
	ULL		nPModels;	/* Maximum number of probability models */
	unsigned	deltaNum;	/* Numerator of delta */
	unsigned	deltaDen;	/* Denominator of delta */
	unsigned	maxCount;	/* Counters /= 2 if one counter >= maxCount */
	unsigned	mode;
	Threshold	threshold;	/* Threshold of context template */
	HashTable	hTable;
	Array		array;
	}
CModel;

typedef struct
	{
	int size;
	MAFImgCoords *position;
	}
CTemplate;

PModel *CreatePModel(unsigned nSymbols);
FloatPModel *CreateFloatPModel(unsigned nSymbols);
void UpdateCModelCounter(CModel *cModel, unsigned pModelIdx, unsigned symbol);
void RemoveCModelCounter(CModel *cModel, unsigned pModelIdx, unsigned char symbol);
CModel *CreateCModel(unsigned maxCtxSize, unsigned nSymbols, unsigned
  nCtxSymbols, unsigned deltaNum, unsigned deltaDen, unsigned maxCount,
  unsigned hSize, unsigned threshold);
void FreeCModel(CModel *);
double FractionOfPModelsUsed(CModel *cModel);
double FractionOfPModelsUsedOnce(CModel *cModel);
void ComputePModel(CModel *cModel, PModel *pModel, unsigned pModelIdx);
void ComputeStaticPModel(MAFImg *mafImg, int row, int col, PModel *pModel);
double PModelSymbolNats(PModel *pModel, unsigned symbol);
void HashingStats(CModel *cModel);
//int PModelUsed(CModel *cModel, unsigned pModelIdx);
int GetPModelIdx(MAFImg *mafImg, int row, int col, CModel *cModel,
  CTemplate *cTemplate);
int GetPModelIdx2(MAFImg *mafImg, int row, int col, CModel *cModel,
  CTemplate *cTemplate, char alm);
int GetPModelIdx3(MAFImg *mafImg, int row, int col);
CTemplate *InitTemplate(int templateId);
CTemplate *InitAncestorTemplate(int leftSize, int rigthSize);
CTemplate *InitAncestorTemplate2(int leftSize, int rigthSize);
void ShowTemplate(CTemplate *cTemplate);


#endif /* CONTEXT_H_INCLUDED */


