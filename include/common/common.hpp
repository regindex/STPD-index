// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  common.hpp: Common definitions and functions file
 */

#ifndef COMMON__HPP_
#define COMMON__HPP_

#include <sys/stat.h>
#include <cassert>

#include <sdsl/int_vector.hpp>

#ifndef M64
	#define M64 0
#endif
#if M64
    typedef uint64_t uint_t;
    typedef int64_t  int_t;
#else
    typedef uint32_t uint_t;
    typedef int32_t  int_t;
#endif
typedef int64_t safe_t;
typedef uint64_t usafe_t;
typedef bool bool_t;
typedef char char_t;
typedef uint8_t uchar_t;

#define DEF_BUFFER_SIZE 5000
#define SIGMA 128
#define SIGMA_DNA 4
#define STORE_SIZE 5

static unsigned char dna_to_code_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 5
};

static unsigned char code_to_dna_table[4] = {'A','C','G','T'};

usafe_t get_5bytes_uint(uchar_t *a, uint_t i=0)
{
  uint_t offset = (i+1)*5-1;
  usafe_t ai = 0;
  for(uint_t j=0;j<5;j++) ai = (ai << 8) | a[offset-j];
    
  return ai;
}

uint8_t bitsize(usafe_t x)
{
    if(x == 0) return 1;
    return 64 - __builtin_clzll(x);
}

#endif 