/*------------------------------------------------------------------------------

Copyright 2009-2012 Luis M. O. Matos (luismatos@ua.pt), All Rights Reserved.

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
#include <limits.h>
#include <time.h>
#include "defs.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"
#include "mem.h"
#include "mafImg.h"
#include "context.h"
#include "common.h"

/*----------------------------------------------------------------------------*/

unsigned hSizes[] = {9999991,19999999,29999999,39999983};

/*----------------------------------------------------------------------------*/

void printHelp()
	{
	fprintf(stderr, "More info: The Ancestral Context Template is a line context configuration that contains the\n");
	fprintf(stderr, "           the more frequent symbol per column. The left and right size define the size of\n");
	fprintf(stderr, "           the template starting from the current column. leftSize and rightSize maximum values\n");
	fprintf(stderr, "           allowed are %d and %d respectably.\n\n", MAX_LEFT_ANCESTRAL_SIZE, MAX_RIGHT_ANCESTRAL_SIZE);

	fprintf(stderr, "Static Column Model: This is a special static model where in each symbol to be encoded, a pModel is\n");
	fprintf(stderr, "created based on the statistics of the current column (Probability Mass Function) of each symbol\n");
	fprintf(stderr, "in case of the first row, we decided that the previous pixel (col-1) should be considered as the most\n");
	fprintf(stderr, "frequent symbol in the current column.\n\n");	

	fprintf(stderr, "Columnwise Model 1: It's model where we will consider the most frequent symbol per column to\n");
	fprintf(stderr, "index a table with z rows and z columns, where z denotes the number of symbols.\n\n");
	
	fprintf(stderr, "Columnwise Model n: It's similar to the previous one. The only diference is that we will consider\n");
	fprintf(stderr, "a table with more rows but with the same number of columns. To get the index we will consider the\n");
	fprintf(stderr, "occurence of each symbol per column. First the symbol that occurs more often until the others that\n");
	fprintf(stderr, "are less frequent.\n");
	}

int main(int argc, char *argv[])

	{
	FILE *inpFp, *txtFp, *fpOut = NULL;
	int n, start, size, srcSize, mafImgRowSize = 0, base, c, c2, cModel, mafPixel, 
	  nCModels = 0, s, bestCModel = 0, row, col, templateId,
	  //nCModels = 0, bestWeightCModel, s, bestCModel = 0, row, col, templateId,
	  deltaNum = DEFAULT_PMODEL_DELTA_NUM, deltaDen = DEFAULT_PMODEL_DELTA_DEN,
          maxCount = DEFAULT_PMODEL_MAX_COUNT, *templateIDs = NULL, leftSize, 
	  rightSize, acModelId = -1, cModelIds[] = {-1,-1,-1};
	//char verbose, fName[256], line[LINE_SIZE], src[32], strand, estimation, acm, 
	char verbose,  *line, src[32], strand, estimation, acm, 
	  alm, scm, cm1, cmn, sizeFlag = 'n', deltaFlag = 'n', thresholdFlag = 'n';
	UChar *mafImgRow = NULL, *mafImgRow2 = NULL;
	MAFImg *mafImg = CreateMAFImg();
	double bestWeight, nats = 0, totalWeight = 0, gamma, *cModelNats, bestNats = 0,
	  *cModelWeight, *cModelProb, bestCModelNats, *cModelGlobalNats, totalNats = 0, 
          bestTotalNats = 0, cpuTimeUsed;
	unsigned pModelIdx = 0, *cModelGlobalUsage, threshold = DEFAULT_THRESHOLD, 
	  *thresholds = NULL, totalBlocks = 0;
	//unsigned int info[MAF_MAX_NROWS][2], symbolsInfo[N_SYMBOLS][2], 
	double info[MAF_MAX_NROWS][2], symbolsInfo[N_SYMBOLS][2];
	unsigned int **usageInfo = NULL;
	unsigned long long bytesOutputed = 0, total[2], totalPixels = 0;
	clock_t tic, tac, start_t;

	CModel **cModels;
	CTemplate **cTemplate;
	PModel **pModels, *mxPModel;
	FloatPModel *floatPModel;
	
	start_t = clock();
	
	total[0] = total[1] = 0;
	verbose = 'n';
	estimation = 'n';	// By default the program creates the binary file
	scm = 'n';
	cm1 = 'n';
	cmn = 'n';
	txtFp = stdout;
	gamma = ((int)(0.94 * GAMMA_K)/(double)GAMMA_K);
	acm = 'n';		// Ancestral Context Template Mode 
	alm = 'n';		// Ancestral Line Mode
	
	if(argc < 2)
		{
		fprintf(stderr, "Usage: SACOe [ -o <outFile>]\n");
		fprintf(stderr, "             [ -v ]\n");
		fprintf(stderr, "             [ -h ]\n");
		fprintf(stderr, "             [ -e ] (estimation only, does not create the binary file)\n");
		fprintf(stderr, "             [ -alm ] (Acenstral Line Mode)\n");
		fprintf(stderr, "             [ -scm ] (Static Column Model)\n");
		fprintf(stderr, "             [ -cm1 <n>/<d> t=<threshold> ] (Columnwise Model 1)\n");
		fprintf(stderr, "             [ -cmn <n>/<d> t=<threshold> ] (Columnwise Model n)\n");
		fprintf(stderr, "             [ -u 0 leftSize-rightSize <n>/<d> t=<threshold> ]" 
		  " (Ancestral Context Template)\n");
		fprintf(stderr, "             [ -u template <n>/<d> t=<threshold> ]\n");
		fprintf(stderr, "             [ -u template <n>/<d> t=<threshold> ]\n");
		fprintf(stderr, "              ...\n");
		fprintf(stderr, "             [ -g | -a <gamma> (def %g) ]\n", gamma);
		fprintf(stderr, "             [ -mc <maxCount> (def %d) ]\n", maxCount);
		fprintf(stderr, "              <mafFile>\n\n");
		return 1;
		}
	
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-h", argv[n]) == 0)
			{
			printHelp();
			return 1;
			}
		
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-e", argv[n]) == 0)
			{
			estimation = 'y';
			break;
			}
	
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-alm", argv[n]) == 0)
			{
			alm = 'y';
			break;
			}
			
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-scm", argv[n]) == 0)
			{
			scm = 'y';
			nCModels++;
			break;
			}
			
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-cm1", argv[n]) == 0)
			{
			cm1 = 'y';
			nCModels++;
			break;
			}

	for(n = 1 ; n < argc ; n++)
		if(strcmp("-cmn", argv[n]) == 0)
			{
			cmn = 'y';
			nCModels++;
			break;
			}
	
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-o", argv[n]) == 0)
			{
			if(estimation == 'y')
				{
				fprintf(stderr, "Attention: estimation mode will be disabled!\n" );
				estimation = 'n';
				}
			if(!(fpOut = fopen(argv[n+1], "w")))
				{
				fprintf(stderr, "Error: can't open file %s\n", argv[n+1]);
				return 1;
				}
			break;
			}
	
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-v", argv[n]) == 0)
			{
			verbose = 'y';
			break;
			}
	
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-g", argv[n]) == 0 || strcmp("-a", argv[n]) == 0)
			{
			gamma = ((int)(atof(argv[n+1]) * GAMMA_K)/(double)GAMMA_K);
			if(gamma >= 1 || gamma < 0)
				{
				fprintf(stderr, "Error: gamma should belong to [0, 1)\n");
				return 1;
				}
			break;
			}

	for(n = 1 ; n < argc ; n++)
		if(strcmp("-mc", argv[n]) == 0)
			{
			maxCount = atoi(argv[n+1]);
			break;
			}
	
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-u", argv[n]) == 0)
			nCModels++;
		
	if(nCModels == 0)
		{
		fprintf(stderr,"Error: at least one model has to be specified\n");
		return 1;
		}

	cModels = (CModel **)Calloc(nCModels, sizeof(CModel *));
	cTemplate = (CTemplate **)Calloc(nCModels, sizeof(CTemplate *));
	pModels = (PModel **)Calloc(nCModels, sizeof(PModel *));
	cModelNats = (double *)Calloc(nCModels, sizeof(double));
	cModelWeight = (double *)Calloc(nCModels, sizeof(double));
	cModelProb = (double *)Calloc(nCModels, sizeof(double));
	cModelGlobalNats = (double *)Calloc(nCModels, sizeof(double));
	cModelGlobalUsage = (unsigned *)Calloc(nCModels, sizeof(unsigned));
	templateIDs = (int *)Calloc(nCModels, sizeof(int));
	thresholds = (unsigned int *)Calloc(nCModels, sizeof(unsigned int));
	usageInfo = (unsigned int **)Calloc(MAF_MAX_NROWS, sizeof(unsigned int *));
	line = (char *)Calloc(LINE_SIZE, sizeof(char));
	
	for(n = 0 ; n < MAF_MAX_NROWS ; n++)
		usageInfo[n] = (unsigned int *)Calloc(nCModels, sizeof(unsigned int));
	
	for(n = 0; n < MAF_MAX_NROWS; n++)
		info[n][0] = info[n][1] = 0.0;
	
	for(n = 0; n < N_SYMBOLS; n++)
		symbolsInfo[n][0] = symbolsInfo[n][1] = 0.0;
		
	cModel = 0;
	for(n = 1 ; n < argc ; n++)
		{
		if(strcmp("-u", argv[n]) == 0)
			{
			templateId = atoi(argv[n+1]);
			
			// Typical template 
			if(templateId != 0)
				{
				cTemplate[cModel] = InitTemplate(templateId);
				if(sscanf(argv[n+2], "%d/%d", &deltaNum, &deltaDen) != 2)
					{
					deltaNum = DEFAULT_PMODEL_DELTA_NUM;
					deltaDen = DEFAULT_PMODEL_DELTA_DEN;
					}
				
				if(argc > n+3 && sscanf(argv[n+2], "t=%u", &threshold) != 1 && 
				  sscanf(argv[n+3], "t=%u", &threshold) != 1)
	                                threshold = DEFAULT_THRESHOLD;
				}
			// Ancestral template
			else
				{
				// Acenstral context mode
				acm = 'y';
				acModelId = cModel;

				// Loop all the parameters for the Ancestral context mode
				for(c=n+2; c < argc-1; c++)
					{
					// Found another template, break the cycle
					if(strcmp(argv[c], "-u") == 0)
						break; // Break the cycle
					else
						{
						// Ancestral template size 
						if(sizeFlag == 'n')
							if(sscanf(argv[c], "%d-%d", &leftSize, &rightSize) == 2)
								sizeFlag = 'y';
						
						if(deltaFlag == 'n')
							if(sscanf(argv[c], "%d/%d", &deltaNum, &deltaDen)  == 2)
								deltaFlag = 'y';
						
						if(thresholdFlag == 'n')
							if(sscanf(argv[n+2], "t=%u", &threshold)  == 1)
								thresholdFlag = 'y';
						}
					}
				
				// Set the default size
				if(sizeFlag == 'n')
					{
					leftSize = DEFAULT_LEFT_ANCESTRAL_SIZE;
					rightSize = DEFAULT_RIGHT_ANCESTRAL_SIZE;
					}

				if(leftSize > MAX_LEFT_ANCESTRAL_SIZE) 
					leftSize = MAX_LEFT_ANCESTRAL_SIZE;

				if(rightSize > MAX_RIGHT_ANCESTRAL_SIZE) 
					rightSize = MAX_RIGHT_ANCESTRAL_SIZE;
									
				if(deltaFlag == 'n')
					{
					deltaNum = DEFAULT_PMODEL_DELTA_NUM;
					deltaDen = DEFAULT_PMODEL_DELTA_DEN;
					}
				
				if(thresholdFlag == 'n')
					threshold = DEFAULT_THRESHOLD;
				
				cTemplate[cModel] = InitAncestorTemplate(leftSize, rightSize);
				}
			
			printf("Using template:\n\n");
			ShowTemplate(cTemplate[cModel]);
			putchar('\n');
			
			cModels[cModel] = CreateCModel(cTemplate[cModel]->size, N_SYMBOLS,
			  N_CTX_SYMBOLS, deltaNum, deltaDen, maxCount, hSizes[0], threshold);

 			if(threshold == 0)
				printf("Creating %"PRIu64" probability models (template size: %d, delta "
                        	  "= %d/%d)\n", (uint64_t)cModels[cModel]->nPModels, 
 	                    	  cModels[cModel]->ctxSize, deltaNum, deltaDen);
			else
				printf("Creating %"PRIu64" probability models (template size: %d, delta "
				  "= %d/%d, threshold = %u)\n", (uint64_t)cModels[cModel]->nPModels, 
				  cModels[cModel]->ctxSize, deltaNum, deltaDen, threshold);

			cModelProb[cModel] = 1.0 / nCModels;
			cModelWeight[cModel] = cModelProb[cModel];
			totalWeight += cModelWeight[cModel];
			pModels[cModel] = CreatePModel(N_SYMBOLS);
			templateIDs[cModel] = templateId;
			thresholds[cModel] = threshold;
			cModel++;
			}
		
		// Static Column Model
		if(strcmp("-scm", argv[n]) == 0)
			{
			// It's not necessary for the Static Column Model
			cModels[cModel] = NULL;
			printf("Creating Static Column Model\n");
			cModelIds[0] = cModel;
			cModelProb[cModel] = 1.0 / nCModels;
			cModelWeight[cModel] = cModelProb[cModel];
			totalWeight += cModelWeight[cModel];
			pModels[cModel] = CreatePModel(N_SYMBOLS);
			cModel++;
			}

		// Columnwise Model 1
		if(strcmp("-cm1", argv[n]) == 0)
			{
			if(sscanf(argv[n+1], "%d/%d", &deltaNum, &deltaDen) != 2)
				{
				deltaNum = DEFAULT_PMODEL_DELTA_NUM;
				deltaDen = DEFAULT_PMODEL_DELTA_DEN;
				}
			
			if(argc > n+3 && sscanf(argv[n+1], "t=%u", &threshold) != 1 && 
			  sscanf(argv[n+2], "t=%u", &threshold) != 1)
				threshold = DEFAULT_THRESHOLD;
			
			// This table is indexed using a single symbol
			cModels[cModel] = CreateCModel(1, N_SYMBOLS, N_CTX_SYMBOLS, 
			  deltaNum, deltaDen, maxCount, hSizes[0], threshold);
				
			if(threshold == 0)
				printf("Creating %"PRIu64" probability models (size: %d, delta "
	                          "= %d/%d) - Columnwise Model 1\n", (uint64_t)cModels[cModel]->nPModels, 
				  cModels[cModel]->ctxSize, deltaNum, deltaDen);
			else
				printf("Creating %"PRIu64" probability models (size: %d, delta "
				  "= %d/%d, threshold = %u) - Columnwise Model 1\n", 
				  (uint64_t)cModels[cModel]->nPModels, cModels[cModel]->ctxSize, deltaNum, 
				  deltaDen, threshold);
				
			cModelIds[1] = cModel;
			cModelProb[cModel] = 1.0 / nCModels;
			cModelWeight[cModel] = cModelProb[cModel];
			totalWeight += cModelWeight[cModel];
			pModels[cModel] = CreatePModel(N_SYMBOLS);
			cModel++;
			}

		// Columnwise Model n
		if(strcmp("-cmn", argv[n]) == 0)
			{
			if(sscanf(argv[n+1], "%d/%d", &deltaNum, &deltaDen) != 2)
				{
				deltaNum = DEFAULT_PMODEL_DELTA_NUM;
				deltaDen = DEFAULT_PMODEL_DELTA_DEN;
				}
			
			if(argc > n+3 && sscanf(argv[n+1], "t=%u", &threshold) != 1 && 
			  sscanf(argv[n+2], "t=%u", &threshold) != 1)
				threshold = DEFAULT_THRESHOLD;
			
			// This table is indexed using a single symbol
			cModels[cModel] = CreateCModel(N_CTX_SYMBOLS, N_SYMBOLS, N_CTX_SYMBOLS, 
			  deltaNum, deltaDen, maxCount, hSizes[0], threshold);
				
			if(threshold == 0)
				printf("Creating %"PRIu64" probability models (size: %d, delta "
	                          "= %d/%d) - Columnwise Model n\n", (uint64_t)cModels[cModel]->nPModels, 
				  cModels[cModel]->ctxSize, deltaNum, deltaDen);
			else
				printf("Creating %"PRIu64" probability models (size: %d, delta "
			          "= %d/%d, threshold = %u) - Columnwise Model n\n", 
				  (uint64_t)cModels[cModel]->nPModels, cModels[cModel]->ctxSize, deltaNum, 
				  deltaDen, threshold);
				
			cModelIds[2] = cModel;
			cModelProb[cModel] = 1.0 / nCModels;
			cModelWeight[cModel] = cModelProb[cModel];
			totalWeight += cModelWeight[cModel];
			pModels[cModel] = CreatePModel(N_SYMBOLS);
			cModel++;
			}
		}
	
	mxPModel = CreatePModel(N_SYMBOLS);
	floatPModel = CreateFloatPModel(N_SYMBOLS);
		
	// If there is a problem in the fopen operation
	if(fpOut == NULL)
		fpOut = fopen("/dev/null", "w");
		
	// Find out if the file can be open for reading
	if(!(inpFp = fopen(argv[argc-1], "r")))
		{
		fprintf(stderr, "Error: unable to open file\n");
		return 1;
		}

	fclose(inpFp);

	/*
	// Open the file for possible decompression with gzip
	#if defined(WIN32) || defined(MSDOS) || defined(__APPLE__) 
		sprintf(fName, "gzip -c -d %s", argv[argc-1]); 
	#else
		sprintf(fName, "zcat -f %s", argv[argc-1]);
	#endif

	if((inpFp = popen(fName, "r")) == NULL)
		{
		fprintf(stderr, "Error: unable to open file\n");
		return 1;
		}
	
	*/

	if((inpFp = fopen(argv[argc - 1], "r")) == NULL)
		{
		fprintf(stderr, "Error: unable to open file\n");
		return 1;
		}
	
	tic = clock();
	
	if(estimation == 'n')
		{
		startoutputtingbits();
		start_encode();

		WriteNBits(maxCount, STORAGE_BITS_PMODEL_MAX_COUNT, fpOut);
		WriteNBits(hSizes[0], STORAGE_BITS_HSIZE, fpOut);
		WriteNBits(nCModels - 1, STORAGE_BITS_N_CMODELS, fpOut);
		WriteNBits(((int)(gamma*GAMMA_K)), STORAGE_BITS_GAMMA, fpOut);
		
		// Ancestral Line Mode
		WriteNBits(alm == 'y' ? 1 : 0, 1, fpOut);
		
		// Column Models mode ID's
		WriteNBits(((scm == 'y' && cModelIds[0] >= 0) ? cModelIds[0]+1:0), 
		  STORAGE_BITS_TEMPLATE_ID, fpOut);
		WriteNBits(((cm1 == 'y' && cModelIds[1] >= 0) ? cModelIds[1]+1:0), 
		  STORAGE_BITS_TEMPLATE_ID, fpOut);
		WriteNBits(((cmn == 'y' && cModelIds[2] >= 0) ? cModelIds[2]+1:0), 
		  STORAGE_BITS_TEMPLATE_ID, fpOut);
		
		for(cModel = 0 ; cModel < nCModels ; cModel++)
			{
			// Columnwise Model 1
			if(cModelIds[1] == cModel || cModelIds[2] == cModel)
				{					
				WriteNBits(cModels[cModel]->deltaNum, STORAGE_BITS_PMODEL_DELTA_NUM, fpOut);
				WriteNBits(cModels[cModel]->deltaDen, STORAGE_BITS_PMODEL_DELTA_DEN, fpOut);
				WriteNBits(thresholds[cModel], STORAGE_BITS_THRESHOLD, fpOut);
				}
			//else	// Write the information for the tipical template
			else if(cModelIds[0] != cModel)
				{
				WriteNBits(templateIDs[cModel], STORAGE_BITS_TEMPLATE_ID, fpOut);
				// Acenstral context? 
				if(acm == 'y' && acModelId == cModel)
					{
					// Storing the left and right size of the ancestral context 
					WriteNBits(leftSize, STORAGE_BITS_LEFT_ANCESTRAL_SIZE,fpOut);
					WriteNBits(rightSize, STORAGE_BITS_RIGHT_ANCESTRAL_SIZE,fpOut);
					//printf("Left Size: %d | Right Size: %d\n", leftSize, rightSize);
					}

				WriteNBits(cModels[cModel]->deltaNum, STORAGE_BITS_PMODEL_DELTA_NUM, fpOut);
				WriteNBits(cModels[cModel]->deltaDen, STORAGE_BITS_PMODEL_DELTA_DEN, fpOut);
				WriteNBits(thresholds[cModel], STORAGE_BITS_THRESHOLD, fpOut);
				}
			}		
		}


	// Read the initial lines until the first "a" line
	do
		{
		if(!fgets(line, LINE_SIZE, inpFp))
			{
			fprintf(stderr, "Error: unexpected end-of-file\n");
			return 1;
			}

		if(verbose == 'y')
			fprintf(txtFp, "%s", line);

		}
	while(line[0] != 'a');
	
	bytesOutputed = 0;
	
	while((c = fgetc(inpFp)) != EOF)
		{
		switch(c)
		{
			case 's': 
				c2 = 's';
				break;

			case '\n':
				c2 = 'e';
				break;

			case 'q':
			case 'i':
			case 'e':
				do
				{
					// Get the remaing line
					if(!fgets(line, LINE_SIZE, inpFp))
						goto end;
					if(verbose == 'y')
						fprintf(txtFp,"%s", line);
					// Get the first character
					c = fgetc(inpFp);
					if(c == EOF) goto end;
				}while(c == 'q' || c == 'e' || c == 'i');
				//}while(line[0] == 'e' || line[0] == 'i' || line[0] == 'q');
				if(c == '\n')
					c2 = 'e';
				else if(c == 's')
					c2 = 's';
				else
					{		
					fprintf(stderr,"Error: unexpected character \"%c\"\n",c);
					return 1;
					}
				break;
			default:
				c2 = 'e';

		}

		switch(c2)
			{
			case 's': break; // Tag of "s" line found
				
			//case '\n': // End of alignment block
			case 'e': // End of alignment block

				if(verbose == 'y')
					printf("Encoding a %d x %d image\n", mafImg->nRows, mafImg->nCols);
				
				if(estimation == 'n')
					{
					// Store the information regarding to each MAF block size
					WriteNBits(mafImg->nRows, STORAGE_BITS_MAF_ROWS, fpOut);
					WriteNBits(mafImg->nCols, STORAGE_BITS_MAF_COLS, fpOut);
					}

				totalBlocks++;
								
				for(row = ((acm == 'y' || alm == 'y' || cm1 == 'y') ? 1 : 0); 
				  row < mafImg->nRows ; row++)
					{
					for(col = 0 ; col < mafImg->nCols ; col++)
						{
						//bestWeightCModel = 0;
						bestWeight = 0;

						for(s = 0 ; s < N_SYMBOLS ; s++)
							floatPModel->freqs[s] = 0;

						bestCModelNats = DBL_MAX;
						mafPixel = GetMAFPixel(mafImg, row, col);

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
							if(acModelId != cModel && cModelIds[0] != cModel && cModelIds[1] != cModel 
							  && cModelIds[2] != cModel)
								{
								pModelIdx = GetPModelIdx2(mafImg, row, col, cModels[cModel], 
								  cTemplate[cModel], ((alm == 'n' && (acm=='y' || 
								  scm == 'y')) ? 'n' : 'y'));
								ComputePModel(cModels[cModel], pModels[cModel], pModelIdx);
								}
							
							
							nats = PModelSymbolNats(pModels[cModel], mafPixel);
							cModelNats[cModel] = nats;
							if(nats < bestCModelNats)
								{
								bestCModelNats = nats;
								bestCModel = cModel;
								}

							if(cModelWeight[cModel] > bestWeight)
								{
								bestWeight = cModelWeight[cModel];
								//bestWeightCModel = cModel;
								}

							for(s = 0 ; s < N_SYMBOLS ; s++)
								floatPModel->freqs[s] += (double) pModels[cModel]->freqs[s] /
								  pModels[cModel]->sum * (cModelWeight[cModel] / totalWeight);

							}
													
						cModelGlobalUsage[bestCModel]++;
						cModelGlobalNats[bestCModel] += bestCModelNats;
						usageInfo[(acm == 'y' || alm == 'y' || cm1 == 'y') ? 
						  mafImg->nRows-1 : mafImg->nRows][bestCModel]++;
					
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
						
						nats = PModelSymbolNats(mxPModel, mafPixel);
						totalNats += nats;
						bestNats += bestCModelNats;
						
						if(estimation == 'n')
							ArithEncodeSymbol(mafPixel, (int *)mxPModel->freqs, 
							  (int) mxPModel->sum, fpOut);
						
						bestTotalNats += bestCModelNats;
												
						// Update models and weights
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
									pModelIdx = GetPModelIdx(mafImg, 0, col, cModels[cModel], 
									  cTemplate[cModel]);

								// Common template
								if(acModelId != cModel && cModelIds[0] != cModel && cModelIds[1] != 
								  cModel && cModelIds[2] != cModel)
									pModelIdx = GetPModelIdx2(mafImg, row, col, cModels[cModel], 
									  cTemplate[cModel], ((alm == 'n' && (acm=='y' || 
									  scm == 'y')) ? 'n' : 'y'));
								}
							
							cModelProb[cModel] = Pow(cModelProb[cModel],gamma) * 
							  (double)pModels[cModel]->freqs[mafPixel]/pModels[cModel]->sum;
	                        			cModelWeight[cModel] = cModelProb[cModel];
	                        			totalWeight += cModelWeight[cModel];
							
							// The current cModel isn't the Static Column Mode
							if(cModelIds[0] != cModel)
								{
								// Store threshold values and remove counter
                            					if(cModels[cModel]->threshold.sizeThreshold != 0 && 
								  totalPixels > cModels[cModel]->threshold.sizeThreshold)
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

						totalPixels++;
												
						// Update estimated ancestor line 
						if(acm == 'y' || alm == 'y' || cm1 == 'y')
							updateAncestorLine(mafImg, row, col);
						
						symbolsInfo[mafPixel][0] += 1.0;
						symbolsInfo[mafPixel][1] += nats;
						
						}
					}
				
				// Update number of pixels
				info[(acm == 'y' || alm == 'y' || cm1 == 'y') ? mafImg->nRows-1 : 
				  mafImg->nRows][0] += ((acm == 'y' || alm == 'y' || cm1 == 'y') ? 
				  mafImg->nRows-1 : mafImg->nRows) * mafImg->nCols;

				// Update number bytes used
				info[(acm == 'y' || alm == 'y' || cm1 == 'y') ? mafImg->nRows-1 : 
				  mafImg->nRows][1] += (unsigned int)(_bytes_output - bytesOutputed);
				bytesOutputed = _bytes_output;
				
				ResetMAFImg(mafImg);
				do
					{
					if(!fgets(line, LINE_SIZE, inpFp))
						goto end;

					if(verbose == 'y')
						fprintf(txtFp, "%s", line);
					}
				while(line[0] != 'a');

				if((c = fgetc(inpFp)) != 's')
					{
					fprintf(stderr, "Error: unexpected character \"%c\"\n", c);
					return 1;
					}

				break;
			default:
				fprintf(stderr, "Error: unexpected character \"%c\"\n", c);
				return 1;

			}

			if(fscanf(inpFp, "%s %d %d %c %d ", src, &start, &size, &strand,
			  &srcSize) != 5)
				{
				fprintf(stderr, "Error: failed to get the \"s\" line\n");
				return 1;
				}

			if(verbose == 'y' && c == 'q')
				printf("Src: %s, Start: %d, Size: %d, Strand: %c, SrcSize: %d\n",
				  src, start, size, strand, srcSize);
		
			while((base = fgetc(inpFp)) != '\n')
				{
				StoreSymbol(&mafImgRow, mafImgRowSize, base);
				mafImgRowSize++;
				}
		
			// Are we in the first line? The Ancestral Context Mode is activated?
			if(mafImg->nRows == 0 && (acm == 'y' || alm == 'y' || cm1 == 'y'))
				{			
				if(!(mafImgRow2 = (UChar *)Calloc((((mafImgRowSize - 1) / 
				  ROW_BLOCK_SIZE) + 1) * ROW_BLOCK_SIZE, sizeof(UChar))))
					{
					fprintf(stderr, "Error: out of memory\n");
					exit(1);
					}

				// Add an ancestral line before adding the sequences aligments
				AddRowToMAFImg(mafImg, mafImgRow2, mafImgRowSize);
				}
		
			AddRowToMAFImg(mafImg, mafImgRow, mafImgRowSize);
			mafImgRow = NULL;
			mafImgRow2 = NULL;
			mafImgRowSize = 0;
			
		}

	end:
	
	if(estimation == 'n')
		{
		// Sending a zero to sinalize the end of the encoding process
		WriteNBits(0, STORAGE_BITS_MAF_ROWS, fpOut);
		finish_encode(fpOut);
		doneoutputtingbits(fpOut);
		}
	
	printf("+--------+--------------+----------------+");
	//printf("| NLines | Total Pixels | Bits per pixel |\n");
	for(cModel = 0 ; cModel < nCModels ; cModel++)
		printf("--------+");
	printf("\n");
	
	printf("| NLines | Total Pixels | Bits per pixel |");
	for(cModel = 0 ; cModel < nCModels ; cModel++)
		{			
		if(cModel == cModelIds[0])
			printf("  -scm  |");

		if(cModel == cModelIds[1])
			printf("  -cm1  |");

		if(cModel == cModelIds[2])
			printf("  -cmn  |");

		if(cModel == acModelId)
			printf("  -acm  |");

		if(cModel != cModelIds[0] && cModel != cModelIds[1] && 
		  cModel != cModelIds[2] && cModel != acModelId)
			printf("   %02d   |", templateIDs[cModel]);
		}

	printf("\n");
	printf("+--------+--------------+----------------+");
	for(cModel = 0 ; cModel < nCModels ; cModel++)
		printf("--------+");
	printf("\n");
	for(n = 0; n < MAF_MAX_NROWS; n++)
		{
		if(info[n][0] > 0.0)
			{
			//printf("|   %02d   | %12u | %10.3f bpp |\n", n, info[n][0], 
			// (double)(info[n][1]*8.0)/info[n][0]);
			printf("|   %02d   | %12u | %10.3f bpp |", n, (unsigned int)info[n][0], 
			  (double)(info[n][1]*8.0)/info[n][0]);
			for(cModel = 0 ; cModel < nCModels ; cModel++)
				printf(" %6.2f |", ((double)usageInfo[n][cModel]/info[n][0])*100.0);
			printf("\n");
			}	
		}	
	
	printf("+--------+--------------+----------------+");
	for(cModel = 0 ; cModel < nCModels ; cModel++)
		printf("--------+");
	printf("\n");
	
	printf("|  Base  |    Count     | Bits per base  |\n");
	printf("+--------+--------------+----------------+\n");
	for(n = 0; n < N_SYMBOLS; n++)
		{
		printf("|   %-2c   | %12u | %10.3f bpp |\n", SymbolToBase(n), (unsigned int)symbolsInfo[n][0], 
		  (double)(symbolsInfo[n][1])/M_LN2/symbolsInfo[n][0]);
		total[0] += symbolsInfo[n][0];
		total[1] += symbolsInfo[n][1];
		}
	
	total[0] -= symbolsInfo[N_SYMBOLS-1][0];
	total[1] -= symbolsInfo[N_SYMBOLS-1][1]; 

	printf("|  ACGT  | %12u | %10.3f bpp |\n", (unsigned)total[0], 
	  (double)(total[1]/M_LN2/total[0]));
	printf("+--------+--------------+----------------+\n");
	
	printf("| Model  |    Count     |    Percent     |\n");
	printf("+--------+--------------+----------------+\n");
	
	for(cModel = 0 ; cModel < nCModels ; cModel++)
		{			
		if(cModel == cModelIds[0])
			printf("|  -scm  | %12u | %12.3f %% | (Static Column Model)\n", 
			  cModelGlobalUsage[cModel], ((double)cModelGlobalUsage[cModel]/
			  totalPixels)*100.0);

		if(cModel == cModelIds[1])
			printf("|  -cm1  | %12u | %12.3f %% | (Columnwise Model 1)\n", 
			  cModelGlobalUsage[cModel], ((double)cModelGlobalUsage[cModel]/
			  totalPixels)*100.0);		

		if(cModel == cModelIds[2])
			printf("|  -cmn  | %12u | %12.3f %% | (Columnwise Model n)\n", 
			  cModelGlobalUsage[cModel], ((double)cModelGlobalUsage[cModel]/
			  totalPixels)*100.0);

		if(cModel == acModelId)
			printf("|  -acm  | %12u | %12.3f %% | (Ancestral Contex Model)\n", 
			  cModelGlobalUsage[cModel], ((double)cModelGlobalUsage[cModel]/
			  totalPixels)*100.0);

		if(cModel != cModelIds[0] && cModel != cModelIds[1] && 
		  cModel != cModelIds[2] && cModel != acModelId)
			printf("|   %02d   | %12u | %12.3f %% |\n", templateIDs[cModel], 
			  cModelGlobalUsage[cModel], ((double)cModelGlobalUsage[cModel]/
			  totalPixels)*100.0);	
		}
	
	printf("+--------+--------------+----------------+\n");
	
	tac = clock();

	printf("Total memory in use: %"PRIu64" bytes\n", (uint64_t)TotalMemory());
	printf("Total blocks: %d\n", totalBlocks);
	printf("Total pixels: %"PRIu64"\n", (uint64_t)totalPixels);
	printf("Total bits: %0.2f ( %0.4f bps ); Lower bound: %0.4f bps\n", 
	  totalNats / M_LN2, totalNats / M_LN2 / totalPixels, 
	  bestTotalNats/M_LN2/totalPixels);
	printf("TOTAL bits: %0.2f (%0.4f bps)\n", _bytes_output*8.0, 
	  _bytes_output * 8.0 /totalPixels);

	cpuTimeUsed = ((double) (tac - tic)) / CLOCKS_PER_SEC;
    	printf("Needed %g seconds for compressing the MAF file.\n", cpuTimeUsed);	
	cpuTimeUsed = ((double) (tac - start_t)) / CLOCKS_PER_SEC;
    	printf("Total cpu time used: %g seconds\n", cpuTimeUsed);
	
	
	fclose(inpFp);
	fclose(txtFp);
	
	return 0;
	}
