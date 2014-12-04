//----------------------------------------------------------------------------|
//    Copyright 2012 University of Aveiro, Portugal, All Rights Reserved.     |
//    This file is part of MAFEnc3,  contact:  ap@ua.pt                       |
//    Description: memory handling functions                                  |
//----------------------------------------------------------------------------|

#define ROW_BLOCK_SIZE 1024

#ifndef MEM_H_INCLUDED
#define MEM_H_INCLUDED

#include <stdio.h>

void *Malloc(size_t size);
void *Calloc(size_t nmemb, size_t size);
void *Realloc(void *ptr, size_t size, size_t additionalSize);
void Free(void *ptr, size_t size);
unsigned long long TotalMemory();

#endif
