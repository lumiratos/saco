/*------------------------------------------------------------------------------

Copyright 2009-2012 Lui M. O. Matos (luismatos@ua.pt), All Rights Reserved.

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
#include <math.h>
#include <float.h>
#include <time.h>
#include "defs.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"
#include "mem.h"
#include "mafImg.h"
#include "context.h"
#include "common.h"

//-----------------------------------------------------------------------------

void writeMAFImg(MAFImg *mafImg, FILE *fp, char am)
	{
	int row, col;

	// Write the image size
	fprintf(fp,"%dx%d\n", (am=='y' ? mafImg->nRows-1:mafImg->nRows), 
	  mafImg->nCols);
	
	// Write the image data
	for(row = (am=='y' ? 1:0) ; row < mafImg->nRows ; row++)
		{
		for(col = 0 ; col < mafImg->nCols ; col++)
			OutputSymbol(GetMAFPixel(mafImg, row, col),fp);

		fprintf(fp,"\n");
		}
	}

//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
	{
	FILE *inpFp = NULL, *outFp = NULL;
	int deltaNum = DEFAULT_PMODEL_DELTA_NUM, deltaDen = DEFAULT_PMODEL_DELTA_DEN,
	  maxCount = DEFAULT_PMODEL_MAX_COUNT, hSize, nCModels, cModel, templateId, 
	  row, col, nRows, nCols, n, s, mafPixel, acModelId=-1, leftSize, rightSize, 
	  cModelIds[] = {-1,-1,-1}, res = -1;
	char verbose = 'n', acm = 'n', alm = 'n', scm = 'n', cm1 = 'n', cmn = 'n';
	double totalWeight = 0, gamma, *cModelWeight, *cModelProb, cpuTimeUsed;
	unsigned pModelIdx = 0, threshold, totalBlocks = 0, totalPixels = 0;
	clock_t tic, tac, start_t;

	CModel **cModels = NULL;
	CTemplate **cTemplate = NULL;
	PModel **pModels = NULL, *mxPModel = NULL;
	FloatPModel *floatPModel = NULL;
	MAFImg *mafImg = NULL;
		
	start_t = clock();
	
	if(argc < 2)
		{
		fprintf(stderr, "Usage: SACOd [ -o outFile]\n");
		fprintf(stderr, "             [ -v ]\n");
		fprintf(stderr, "              <codeFile>\n");
		return 1;
		}

	for(n = 1 ; n < argc ; n++)
		if(strcmp("-v", argv[n]) == 0)
			{
			verbose = 'y';
			break;
			}

	for(n = 1 ; n < argc ; n++)
		if(strcmp("-o", argv[n]) == 0)
			{
			if(!(outFp = fopen(argv[n+1], "w")))
				{
				fprintf(stderr, "Error: unable to open file %s\n", argv[n+1]);
				return 1;
				}

			break;
			}	

	if(!(inpFp = fopen(argv[argc-1], "r")))
		{
		fprintf(stderr, "Error: unable to open file %s\n", argv[argc-1]);
		return 1;
		}
	
	tic = clock();
	
	startinputtingbits();
	start_decode(inpFp);

	// Read side information
	maxCount = ReadNBits(STORAGE_BITS_PMODEL_MAX_COUNT, inpFp);		
	hSize = ReadNBits(STORAGE_BITS_HSIZE, inpFp);	
	nCModels = ReadNBits(STORAGE_BITS_N_CMODELS, inpFp) + 1;
	gamma = ReadNBits(STORAGE_BITS_GAMMA, inpFp) / (double)GAMMA_K;		
	alm = (ReadNBits(1, inpFp) == 1 ? 'y' :'n');
	
	// Column Models mode ID's - Static Column Model
	res = ReadNBits(STORAGE_BITS_TEMPLATE_ID, inpFp);
	if(res > 0) 
		{
		cModelIds[0] = res - 1;
		scm = 'y';
		}

	// Columnwise Model 1	
	res = ReadNBits(STORAGE_BITS_TEMPLATE_ID, inpFp);
	if(res > 0) 
		{
		cModelIds[1] = res - 1;
		cm1 = 'y';
		}

	// Columnwise Model n		
	res = ReadNBits(STORAGE_BITS_TEMPLATE_ID, inpFp);
	if(res > 0) 
		{
		cModelIds[2] = res - 1;
		cmn = 'y';
		}
	
	// Allocate memory for models and templates
	cModels = (CModel **)Calloc(nCModels, sizeof(CModel *));
	cTemplate = (CTemplate **)Calloc(nCModels, sizeof(CTemplate *));
	pModels = (PModel **)Calloc(nCModels, sizeof(PModel *));
	cModelWeight = (double *)Calloc(nCModels, sizeof(double));
	cModelProb = (double *)Calloc(nCModels, sizeof(double));
	
	// Create the context templates and CModels
	for(cModel = 0 ; cModel < nCModels ; cModel++)
		{
		// For the Columnwise Models (not the static one)
		if(cModelIds[1] == cModel || cModelIds[2] == cModel)
			{
			deltaNum = ReadNBits(STORAGE_BITS_PMODEL_DELTA_NUM, inpFp);
			deltaDen = ReadNBits(STORAGE_BITS_PMODEL_DELTA_DEN, inpFp);
			threshold = ReadNBits(STORAGE_BITS_THRESHOLD, inpFp);
				
			// Columnwise Model 1
			if(cModelIds[1] == cModel)
				{
				// This table is indexed using a single symbol
				cModels[cModel] = CreateCModel(1, N_SYMBOLS, N_CTX_SYMBOLS, deltaNum, 
				  deltaDen, maxCount, hSize, threshold);

				if(threshold == 0)
					printf("Creating %"PRIu64" probability models (size: %d, delta = %d/%d) "
					  "- Columnwise Model 1\n", (uint64_t)cModels[cModel]->nPModels, 
					  cModels[cModel]->ctxSize, deltaNum, deltaDen);
				else
					printf("Creating %"PRIu64" probability models (size: %d, delta = %d/%d, "
					  "threshold = %u) - Columnwise Model 1\n", (uint64_t)cModels[cModel]->nPModels, 
					  cModels[cModel]->ctxSize, deltaNum, deltaDen, threshold);
				}
			else // Columnwise Model n
				{
				cModels[cModel] = CreateCModel(N_CTX_SYMBOLS, N_SYMBOLS, N_CTX_SYMBOLS, 
				  deltaNum, deltaDen, maxCount, hSize, threshold);
				if(threshold == 0)
					printf("Creating %"PRIu64" probability models (size: %d, delta "
					  "= %d/%d) - Columnwise Model n\n", (uint64_t)cModels[cModel]->nPModels, 
					  cModels[cModel]->ctxSize, deltaNum, deltaDen);
				else
					printf("Creating %"PRIu64" probability models (size: %d, delta "
					  "= %d/%d, threshold = %u) - Columnwise Model n\n", 
					  (uint64_t)cModels[cModel]->nPModels, cModels[cModel]->ctxSize, 
					  deltaNum, deltaDen, threshold);
				}
			
			cModelProb[cModel] = 1.0 / nCModels;
			cModelWeight[cModel] = cModelProb[cModel];
			totalWeight += cModelWeight[cModel];
			pModels[cModel] = CreatePModel(N_SYMBOLS);	
			}
		else if(cModelIds[0] != cModel)	// For the other models...
			{
			// Read template ID 
			templateId = ReadNBits(STORAGE_BITS_TEMPLATE_ID, inpFp);
			//printf("Template ID = %d\n", templateId);
			// Typical template 
			if(templateId != 0)
				cTemplate[cModel] = InitTemplate(templateId);
			else
				{
				// Acenstral context mode
				acm = 'y';
				acModelId = cModel;
				// Read the template Ancestral size
				leftSize = ReadNBits(STORAGE_BITS_LEFT_ANCESTRAL_SIZE, inpFp);
				rightSize = ReadNBits(STORAGE_BITS_RIGHT_ANCESTRAL_SIZE, inpFp);				
				//printf("left side: %d | right side: %d\n", leftSize, rightSize);
				cTemplate[cModel] = InitAncestorTemplate(leftSize, rightSize);
				}
			
			printf("Using template:\n\n");
			ShowTemplate(cTemplate[cModel]);
			putchar('\n');
			
			deltaNum = ReadNBits(STORAGE_BITS_PMODEL_DELTA_NUM, inpFp);
			deltaDen = ReadNBits(STORAGE_BITS_PMODEL_DELTA_DEN, inpFp);
			threshold = ReadNBits(STORAGE_BITS_THRESHOLD, inpFp);
			cModels[cModel] = CreateCModel(cTemplate[cModel]->size, N_SYMBOLS, 
			  N_SYMBOLS, deltaNum, deltaDen, maxCount, hSize, threshold);
			
			if(threshold == 0)
				printf("Creating %"PRIu64" probability models (template size: %d, delta "
				  "= %d/%d)\n", (uint64_t)cModels[cModel]->nPModels, cModels[cModel]->ctxSize, 
				  deltaNum, deltaDen);
			else
				printf("Creating %"PRIu64" probability models (template size: %d, delta "
				  "= %d/%d, threshold = %u)\n", (uint64_t)cModels[cModel]->nPModels, 
				  cModels[cModel]->ctxSize, deltaNum, deltaDen, threshold);
		
			cModelProb[cModel] = 1.0 / nCModels;
			cModelWeight[cModel] = cModelProb[cModel];
			totalWeight += cModelWeight[cModel];
			pModels[cModel] = CreatePModel(N_SYMBOLS);
			}
		else if(cModelIds[0] == cModel) // Static Column Mode
			{
			printf("Creating Static Column Model...\n");
			// It's not necessary
			cModels[cModel] = NULL;
			cModelProb[cModel] = 1.0 / nCModels;
			cModelWeight[cModel] = cModelProb[cModel];
			totalWeight += cModelWeight[cModel];
			pModels[cModel] = CreatePModel(N_SYMBOLS);						
			}
		}
	
	mxPModel = CreatePModel(N_SYMBOLS);
	floatPModel = CreateFloatPModel(N_SYMBOLS);
		
	nRows = -1;
	// Read all sequences until the last one
	while(nRows != 0)
		{
		// Read the information regarding to the number of rows and columns
		nRows = ReadNBits(STORAGE_BITS_MAF_ROWS, inpFp);
	
		// End of the compressed file 
		if(nRows == 0) 
			break;
		else 
			nCols = ReadNBits(STORAGE_BITS_MAF_COLS, inpFp);
		
		// Create MAF Image with nRows rows and
		mafImg = CreateMAFImgWith(nRows, nCols);
		
		totalBlocks++;
		
		if(verbose == 'y')
			printf("Decoding MAF block %u with %d rows and %d "
			  "cols...\n",totalBlocks, nRows, nCols);
			
		/* Read the compressed file MAF block by MAF block */
		for(row = ((acm == 'y' || alm == 'y' || cm1 == 'y') ? 1 : 0); 
		  row < mafImg->nRows ; row++)
			{
			for(col = 0 ; col < nCols ; col++)
				{
				for(s = 0 ; s < N_SYMBOLS ; s++)
					floatPModel->freqs[s] = 0;
				
				for(cModel = 0 ; cModel < nCModels ; cModel++)
					{
					// Static Column Model
					if(scm == 'y' && cModelIds[0] == cModel)
						// Count the number of symbols in the current column
						ComputeStaticPModel(mafImg, row, col, pModels[cModel]);
						
					// Columnwise Model 1
					if(cm1 == 'y' && cModelIds[1] == cModel)
						{
						// The line zero, will contain the most frequent symbol per column
						pModelIdx = GetMAFPixel(mafImg, 0, col);
						ComputePModel(cModels[cModel], pModels[cModel], pModelIdx);
						}
													
					// Columnwise Model n
					if(cmn == 'y' && cModelIds[2] == cModel)
						{
						pModelIdx = GetPModelIdx3(mafImg, row, col);
						ComputePModel(cModels[cModel], pModels[cModel], pModelIdx);
						}
					
					// Get IDX for the ancestral line - Line 0
					if(acm == 'y' && acModelId == cModel)
						{
						pModelIdx = GetPModelIdx(mafImg, 0, col, cModels[cModel], 
						  cTemplate[cModel]);								
						ComputePModel(cModels[cModel], pModels[cModel], pModelIdx);
						}
						
					// For common template	
					if(acModelId != cModel && cModelIds[0] != cModel && cModelIds[1] != 
					  cModel && cModelIds[2] != cModel)
						{
						pModelIdx = GetPModelIdx2(mafImg, row, col, cModels[cModel], 
						  cTemplate[cModel], ((alm == 'n' && (acm=='y' || 
						  scm == 'y')) ? 'n' : 'y'));
						ComputePModel(cModels[cModel], pModels[cModel], pModelIdx);
						}
						
					for(s = 0 ; s < N_SYMBOLS ; s++)
						floatPModel->freqs[s] += (double)pModels[cModel]->freqs[s] / 
						  pModels[cModel]->sum * (cModelWeight[cModel] / totalWeight);
					}
					
				// Transform floating probabilities back to integers.
				// This is needed by the arithmetic encoder. It was left
				// here just for maintaining the parallel with a real
				// encoder.
				mxPModel->sum = 0;
				for(s = 0 ; s < N_SYMBOLS ; s++)
					{
					mxPModel->freqs[s] = (unsigned) (floatPModel->freqs[s] * maxCount);
					if(!mxPModel->freqs[s])
						mxPModel->freqs[s]++;
					mxPModel->sum += mxPModel->freqs[s];
					}
				
				// Decode MAF Pixel
				mafPixel = ArithDecodeSymbol(N_SYMBOLS, (int *)mxPModel->freqs, 
				  (int)mxPModel->sum, inpFp);
				
				/* Update models and weights */
				totalWeight = 0;
				for(cModel = 0 ; cModel < nCModels ; cModel++)
					{
					// The current cModel isn't the Static Column Mode
					if(cModelIds[0] != cModel)
						{
						// Columnwise Model 1
						if(cm1 == 'y' && cModelIds[1] == cModel)
							// The line zero, will contain the most frequent symbol per column
							pModelIdx = GetMAFPixel(mafImg, 0, col);
								
						// Columnwise Model n
						if(cmn == 'y' && cModelIds[2] == cModel)
							pModelIdx = GetPModelIdx3(mafImg, row, col);
						
						// Get IDX for the ancestral line - Line 0
						if(acm == 'y' && acModelId == cModel)
							pModelIdx = GetPModelIdx(mafImg, 0, col, cModels[cModel], cTemplate[cModel]);

						// Common template
						if(acModelId != cModel && cModelIds[0] != cModel && cModelIds[1] != 	
						  cModel && cModelIds[2] != cModel)
							pModelIdx = GetPModelIdx2(mafImg, row, col, cModels[cModel], 
							  cTemplate[cModel], ((alm == 'n' && (acm=='y' || 
							  scm == 'y')) ? 'n' : 'y'));
						}

					cModelProb[cModel] = Pow(cModelProb[cModel],gamma) * 
					  (double)pModels[cModel]->freqs[mafPixel] / pModels[cModel]->sum;
					cModelWeight[cModel] = cModelProb[cModel];
                			totalWeight += cModelWeight[cModel];
											
					// The current cModel isn't the Static Column Mode
					if(cModelIds[0] != cModel)
						{
						// Store threshold values and remove counter
						if(cModels[cModel]->threshold.sizeThreshold != 0 && totalPixels > 
						  cModels[cModel]->threshold.sizeThreshold)
							RemoveCModelCounter(cModels[cModel], cModels[cModel]->
					  		  threshold.idxOcc[cModels[cModel]->threshold.indexThreshold], 
					  		  cModels[cModel]->threshold.symbols[cModels[cModel]->
							  threshold.indexThreshold]);
                                        
						// Update and Move forward the index of the threshold (on the loop)
						if(cModels[cModel]->threshold.sizeThreshold != 0)
							{
            						cModels[cModel]->threshold.idxOcc[cModels[cModel]->
							  threshold.indexThreshold] = pModelIdx;
							cModels[cModel]->threshold.symbols[cModels[cModel]->
							  threshold.indexThreshold] = mafPixel;

							if(cModels[cModel]->threshold.indexThreshold == 
							  cModels[cModel]->threshold.sizeThreshold - 1)
								cModels[cModel]->threshold.indexThreshold = 0;
							else
								cModels[cModel]->threshold.indexThreshold++;
							}

						UpdateCModelCounter(cModels[cModel], pModelIdx, mafPixel);
						}
					}
				
				SetMAFPixel(mafImg, row, col, mafPixel);
				totalPixels++;	
							
				// Update estimated ancestor line 
				if(acm == 'y' || alm == 'y' || cm1 == 'y')
					updateAncestorLine(mafImg, row, col);
				}
			}

		writeMAFImg(mafImg, outFp, ((acm == 'y' || alm == 'y' || 
		  cm1 == 'y') ? 'y' : 'n'));
		FreeMAFImg(mafImg);
		mafImg = NULL;
		}
	
	// The file is all decoded
	finish_decode();
	doneinputtingbits();

	printf("Total blocks: %d decoded\n", totalBlocks);
	
	tac = clock();
	
	cpuTimeUsed = ((double) (tac - tic)) / CLOCKS_PER_SEC;
    	printf("Needed %g seconds for decompressing the compressede file.\n", 
	  cpuTimeUsed);	
	cpuTimeUsed = ((double) (tac - start_t)) / CLOCKS_PER_SEC;
    	printf("Total cpu time used: %g seconds\n", cpuTimeUsed);

	fclose(inpFp);
	fclose(outFp);

	return 0;		
	}
