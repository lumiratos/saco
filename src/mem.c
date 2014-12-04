//----------------------------------------------------------------------------|
//    Copyright 2012 University of Aveiro, Portugal, All Rights Reserved.     |
//    This file is part of MAFEnc3,  contact:  ap@ua.pt                       |
//    Description: memory handling functions                                  |
//----------------------------------------------------------------------------|

#include "mem.h"
#include "defs.h"
#include <stdlib.h>


//-----------------------------------------------------------------------------

static uint64_t totalMemory = 0;
//static unsigned long long totalMemory = 0;

//-----------------------------------------------------------------------------

void *Malloc(size_t size)
        {
        void *pointer = malloc(size);

        if(pointer == NULL)
                {
                fprintf(stderr, "Error allocating %"PRIu64" bytes.\n", (uint64_t)size);
                //fprintf(stderr, "Error allocating %zu bytes.\n", size);
                exit(1);
                }

        totalMemory += size;
		//printf("Allocated %zu bytes of memory (Malloc)\n", size);
        return pointer;
        }

//----------------------------------------------------------------------------

void *Calloc(size_t nmemb, size_t size)
        {
        void *pointer = calloc(nmemb, size);

        if(pointer == NULL)
                {
                fprintf(stderr, "Error allocating %"PRIu64" bytes.\n", (uint64_t)size);
                //fprintf(stderr, "Error allocating %zu bytes.\n", size);
                exit(1);
                }

        totalMemory += nmemb * size;
		//printf("Allocated %zu bytes of memory (Calloc)\n", size*nmemb);
        return pointer;
        }

//----------------------------------------------------------------------------

void *Realloc(void* ptr, size_t size, size_t additionalSize)
        {
        void *pointer = realloc(ptr, size);

        if(pointer == NULL)
                {
		
                fprintf(stderr, "Error allocating %"PRIu64" bytes.\n", (uint64_t)size);
                //fprintf(stderr, "Error allocating %zu bytes.\n", size);
                exit(1);
                }

        totalMemory += additionalSize;
		//printf("Reallocated %zu bytes of memory (Realloc)\n", size);
        return pointer;
        }

//----------------------------------------------------------------------------

void Free(void *ptr, size_t size)
	{
	totalMemory -= size;
	free(ptr);
	//printf("Freed %zu bytes of memory (Free)\n", size);
	}

//----------------------------------------------------------------------------

unsigned long long TotalMemory()
        {
        return totalMemory;
        }

//----------------------------------------------------------------------------

