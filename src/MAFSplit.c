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
#include "defs.h"
#include "mem.h"
#include "mafImg.h"
#include "common.h"


int main(int argc, char *argv[])

	{
	FILE *inpFp = NULL, *infoFp = NULL, *msaFp = NULL, *MSAFp = NULL, *MAFFp = NULL; 
	int c2 = 0, c = 0, base = 0, n;
	//int c2 = 0, c = 0, base = 0, mafImgRowSize = 0, n;
	uint64_t start = 0, size = 0, srcSize = 0;
	//unsigned long long start = 0, size = 0, srcSize = 0;
	//char fName[256], line[LINE_SIZE], src[32], strand = ' ';
	char fName[256], *line, src[32], strand = ' ';
	unsigned int totalBlocks = 0;
	//UChar *mafImgRow = NULL;
	MAFImg *mafImg = CreateMAFImg();
	
	if(argc < 3)
		{
		fprintf(stderr, "Usage: MAFSplit [original MAF file]\n\n");
		fprintf(stderr, "                [ -I -i -info file.info] (file with  all header information)\n");
		fprintf(stderr, "                [ -m -msa file.msa] (file with only the original bases and gaps)\n");
		fprintf(stderr, "                [ -M -MSA file.MSA] (file with only ACGT and gaps '-')\n");
		fprintf(stderr, "                [ -MAF file.MAF] (Original file with ACGT and gaps '-')\n\n");
		fprintf(stderr, "Note: The original MAF file must be a file in .maf format.\n");
		fprintf(stderr, "      The program will create an info file (MAF.info) with\n");
		fprintf(stderr, "      all the information regarding to the input MAF file.\n\n");
		fprintf(stderr, "MAF format: https://cgwb.nci.nih.gov/FAQ/FAQformat.html#format5\n");
		return 1;
		}
	
	for(n = 1 ; n < argc ; n++)
		if((strcmp("-I", argv[n]) == 0) || (strcmp("-i", argv[n]) == 0) || (strcmp("-info", argv[n]) == 0))
			{
			if(!(infoFp = fopen(argv[n+1], "w")))
				{
				fprintf(stderr, "Error: unable to open file %s for writing!\n", argv[n+1]);
				return 1;
				}
			break;	
			}
	
	for(n = 1 ; n < argc ; n++)
		if((strcmp("-m", argv[n]) == 0) || (strcmp("-msa", argv[n]) == 0))
			{
			if(!(msaFp = fopen(argv[n+1], "w")))
				{
				fprintf(stderr, "Error: unable to open file %s for writing!\n", argv[n+1]);
				return 1;
				}
			break;	
			}
			
	for(n = 1 ; n < argc ; n++)
		if((strcmp("-M", argv[n]) == 0) || (strcmp("-MSA", argv[n]) == 0))
			{
			if(!(MSAFp = fopen(argv[n+1], "w")))
				{
				fprintf(stderr, "Error: unable to open file %s for writing!\n", argv[n+1]);
				return 1;
				}
			break;	
			}
	for(n = 1 ; n < argc ; n++)
		if(strcmp("-MAF", argv[n]) == 0)
			{
			if(!(MAFFp = fopen(argv[n+1], "w")))
				{
				fprintf(stderr, "Error: unable to open file %s for writing!\n", argv[n+1]);
				return 1;
				}
			break;	
			}
	
	if(infoFp == NULL)
		infoFp = fopen("/dev/null", "w");
	if(msaFp == NULL)
		msaFp = fopen("/dev/null", "w");
	if(MSAFp == NULL)
		MSAFp = fopen("/dev/null", "w");
	if(MAFFp == NULL)
		MAFFp = fopen("/dev/null", "w");
	
	line = (char *)Calloc(LINE_SIZE, sizeof(char));	

	/* Find out if the file can be open for reading */
	if(!(inpFp = fopen(argv[argc-1], "r")))
		{
		fprintf(stderr, "Error: unable to open file %s for reading!\n", argv[argc-1]);
		return 1;
		}

	fclose(inpFp);

	// Open the file for possible decompression with gzip
	#if defined(WIN32) || defined(MSDOS) || defined(__APPLE__) 
		sprintf(fName, "gzip -c -d %s", argv[argc-1]); 
	#else
		sprintf(fName, "zcat -f %s", argv[argc-1]);
	#endif

	if((inpFp = popen(fName, "r")) == NULL)
		{
		fprintf(stderr, "Error: unable to open file %s\n", fName);
		return 1;
		}
	
	// Read the initial lines until the first "a" line
	do
		{
		if(!fgets(line, LINE_SIZE, inpFp))
			{
			fprintf(stderr, "Error: unexpected end-of-file\n");
			return 1;
			}
			// Header lines until the first score line
			fprintf(infoFp, "%s", line);
			fprintf(MAFFp, "%s", line);
		}
	while(line[0] != 'a');

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
						fprintf(infoFp, "%s", line);
						fprintf(MAFFp, "%s", line);

						// Get the first character
						c = fgetc(inpFp);
						fprintf(infoFp, "%c", c);
						fprintf(MAFFp, "%c", c);
						if(c == EOF) goto end;
					}while(c == 'q' || c == 'e' || c == 'i');

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
				totalBlocks++;
								
				// Free memory
				ResetMAFImg(mafImg);
				
				do
				{
					if(!fgets(line, LINE_SIZE, inpFp))
						goto end;
					fprintf(infoFp, "%s", line);
					fprintf(MAFFp, "%s", line);
				} while(line[0] != 'a');
				
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

		if(fscanf(inpFp, "%s %"SCNu64" %"SCNu64" %c %"SCNu64" ", src, &start, &size, &strand, &srcSize) != 5)
		//if(fscanf(inpFp, "%s %llu %llu %c %llu ", src, &start, &size, &strand, &srcSize) != 5)
			{
			fprintf(stderr, "Error: failed to get the \"s\" line\n");
			return 1;
			}
		fprintf(infoFp, "s %-31s %10"PRIu64" %8"PRIu64" %1c %10"PRIu64"\n", src, (uint64_t)start, (uint64_t)size,strand, (uint64_t)srcSize);
		fprintf(MAFFp, "s %-31s %10"PRIu64" %8"PRIu64" %1c %10"PRIu64"\n", src, (uint64_t)start, (uint64_t)size, strand, (uint64_t)srcSize);
		//fprintf(infoFp, "s %-31s %10llu %8llu %1c %10llu\n", src, start, size, strand, srcSize);
		//fprintf(MAFFp, "s %-31s %10llu %8llu %1c %10llu ", src, start, size, strand, srcSize);
	
			
		while((base = fgetc(inpFp)) != '\n')
			{
			//StoreSymbol(&mafImgRow, mafImgRowSize, base);
			//mafImgRowSize++;
			fprintf(msaFp, "%c", base);
			fprintf(MSAFp, "%c", BaseTransform(base));
			fprintf(MAFFp, "%c", BaseTransform(base));
			}
		fprintf(msaFp, "\n");
		fprintf(MSAFp, "\n");
		fprintf(MAFFp, "\n");

		//AddRowToMAFImg(mafImg, mafImgRow, mafImgRowSize);
		//mafImgRow = NULL;
		//mafImgRowSize = 0;
		}

	end:
			
	fclose(inpFp);
	fclose(infoFp);
	fclose(msaFp);
	fclose(MSAFp);
	fclose(MAFFp);

	Free(line, LINE_SIZE*sizeof(char));
	
	printf("Total blocks: %u\n", totalBlocks);
	printf("Total memory in use: %"PRIu64" bytes.\n", (uint64_t)TotalMemory());
	
	return 0;
	}
