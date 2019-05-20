/*
 * Copyright 2019, José-Manuel Herruzo <jmherruzo@uma.es>,
 *                 Jesús Alastruey-Benedé <jalastru@unizar.es>,
 *                 Pablo Ibáñez-Marín <imarin@unizar.es>
 *
 * This file is part of the bvSFM sequence alignment package.
 *
 * bvSFM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * bvSFM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bvSFM. If not, see <http://www.gnu.org/licenses/>.
 *
 * If you publish any work that uses this software, please cite the following paper:
 *
 * J.M. Herruzo, S. González-Navarro, P. Ibáñez, V. Viñals, J. Alastruey-Benedé, and Óscar Plata.
 * Accelerating Sequence Alignments Based on FM-Index Using the Intel KNL Processor.
 * IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB 2019).
 * DOI: 10.1109/TCBB.2018.2884701 
 * 
 * @article{herruzo2019TCBB,
 *  author    = {José Manuel Herruzo, Sonia González-Navarro, Pablo Ibáñez, Víctor Viñals, Jesús Alastruey-Benedé, and Óscar Plata},
 *  journal = {IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB 2019)},
 *  title     = {Accelerating Sequence Alignments Based on FM-Index Using the Intel KNL Processor},
 *  year      = {2019},
 *  doi       = {10.1109/TCBB.2018.2884701}
 * }
 *
 */

#ifndef _SFM_H_
#define _SFM_H_

#include <immintrin.h>
#include <stdint.h>

typedef unsigned int uint;

#ifdef __GNUC__
#define __forceinline inline
//#define __forceinline __attribute__((always_inline))
#define _popcnt64 __builtin_popcountll
#endif

// k-steps
#define KSTEPS 2

// sampling factor
#define D_VAL    64

// bits to encode a symbol (A,C,G,T)
#define BITS_PER_SYMBOL    2
// #define BITS_PER_K2_SYMBOL 4

// number of symbols
#define SYMBOLS     4
#define K2_SYMBOLS 16

typedef struct SFM_entry {
  uint64_t counter;    // 8 bytes, uint64_t instead of uint32_t (padding)
  // uint32_t padding;
  uint64_t data;       // 8 bytes
} SFM_entry_t;

typedef struct LUT_entry {
  uint32_t start;
  uint32_t end;
} LUT_entry_t;

typedef struct SFM_Index {
  uint64_t len;     // BWT lenght
  char * start;     // first 500-char of raw data
  char * alphabet;  // unique symbols: ACGT
  char last_char;
  uint64_t * end_char_pos;
  uint64_t * C;
  uint64_t n_entries;
  SFM_entry_t * entries;
  uint8_t * encoding_table;   // 1-char step encoding table
  uint8_t * encoding_table2;  // 2-char step encoding table
  LUT_entry_t * LUT[KSTEPS];  // 5 and 6 characters LUTs
} SFM_t;

void init_C(uint64_t C[KSTEPS][SYMBOLS]);

void dump_C(uint64_t C[KSTEPS][SYMBOLS]);

/* nrows: KSTEPS */
void dump_BWT(char **matrix, uint nrows, uint ncols);
void dump_encoded_BWT(char **matrix, uint nrows, uint ncols, char *codes, uint64_t* end);

/* convert from encoded symbols to array of chars */
/* For example:
 * 10 00 11 11 00 01 00
 *  G  A  T  T  A  C  A */
void decode_symbols(uint idx, char* seq_out, char* alphabet, uint bits_symbol, uint len);

void build_C(uint64_t *C, char *data, uint64_t data_len);

void dump_array(char *array, uint64_t len);

void mask_init(uint64_t *mask);

int generate_SFM(SFM_t *fmi, uint8_t** bwt, uint64_t len, char* alphabet, size_t alignment, uint64_t* end_char);

/**
  @param file Char array containing the filename
  @param fmi FMIndex which will be written to the file
  @return 0 if no error appeared.
*/
int write_SFM(const char* file, SFM_t *fmi);

int generate_SFM_LUT(SFM_t * fmi, uint n_chars, LUT_entry_t ** LUT);

int dump_SFM(SFM_t *fmi);

int count_SFM(SFM_t *fmi, const char* orig_seq, uint len, uint64_t * start, uint64_t* end);

/**
  @param file Char array containing the filename
  @param fmi FMIndex to load
  @return 0 if no error appeared.
*/
int load_SFM(const char* file, SFM_t * fmi);

int generate_SFM_encoding_table(SFM_t *fmi, uint8_t** out);

int generate_SFM_encoding_table2(SFM_t *fmi, uint8_t** out);

void free_SFM(SFM_t * fmi);

#endif
