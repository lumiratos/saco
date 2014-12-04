//----------------------------------------------------------------------------|
//    Copyright 2012 University of Aveiro, Portugal, All Rights Reserved.     |
//    This file is part of MAFEnc,  contacts:  luismatos@ua.pt or ap@ua.pt    |
//    Description: defining macros file                                       |
//----------------------------------------------------------------------------|

#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define ROW_BLOCK_SIZE 1024
#define N_SYMBOLS 5
#define N_ALL_SYMBOLS 256
#define N_CTX_SYMBOLS 5
//#define LINE_SIZE 256
//#define LINE_SIZE 16777216			// 2^24
#define LINE_SIZE 16777216			//
#define GAMMA_K	65536
#define DEFAULT_PMODEL_DELTA_NUM	1
#define DEFAULT_PMODEL_DELTA_DEN	1
#define DEFAULT_THRESHOLD		0
#define DEFAULT_LEFT_ANCESTRAL_SIZE 4
#define DEFAULT_RIGHT_ANCESTRAL_SIZE 3
#define MAX_LEFT_ANCESTRAL_SIZE 8
#define MAX_RIGHT_ANCESTRAL_SIZE 8
#define MAF_MAX_NROWS	110	
//#define MAF_MAX_NCOLS	16777216 	// 2^24
#define MAF_MAX_NCOLS	16777216 	//

#endif
