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

#ifndef _TYPES_H_
#define _TYPES_H_

#include <immintrin.h>
#include <stdint.h>

typedef unsigned int uint;

#ifdef __GNUC__
#define __forceinline inline
//#define __forceinline __attribute__((always_inline))
#define _popcnt64 __builtin_popcountll
#endif

struct block {
  uint64_t counter;
  //uint32_t padding;
  uint64_t bitmap;
};


struct FMIndex {
  uint64_t len;
  char * start;

  uint char_count; // number of unique symbols (4)
  uint char_bits;  // bits required to encode the symbols (2)
  char * characters;

  char last_char;
  uint64_t * end_char;
  uint64_t * C;

  uint64_t n_entries;
  uint char_per_entry;

  struct block * entries;

  // 1-char step encoding table
  uint8_t * encoding_table;

  // 2-char step encoding table
  uint8_t * encoding_table2;

  uint *** LUT; // 5 and 6 characters LUTs
};

#endif
