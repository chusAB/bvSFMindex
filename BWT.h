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

#ifndef _BWT_H_
#define _BWT_H_

#include "types.h"


/**
  @param in Char array containing all the data
  @param out Char array containing the distinct elements in In. Allocated
  inside the function
  @result Number of unique elements. Negative if some error occurred.
*/
int get_unique_elements(char * in, char ** out, uint n);

/**
  @param in Char array containing the original text
  @param out Char array containing the BWT transform
  @param n Number of characters
  @param steps Number of steps proccessed in each iteration (see N-Steps FM-Index)
  @result Position of the final character of the string ($). Negative number if
  error.
*/
int get_bwt(char** in, char*** out, uint64_t n, uint steps, uint64_t** end);

/**
  @param bwt Char array containing the BWT. Will be modified inside.
  @param codes Char array containing the character substitution. The character
  at position i will be switched by the i number.
  @param n_bwt Number of characters in BWT
  @param n_characters Number of unique characters
  @result 0 if no error occurred
*/
int encode_bwt(char** bwt, char* codes, uint n_bwt, uint n_characters, uint steps, uint64_t* end);

/**
  @param bwt Char array containing the BWT.
  @param n Number of characters in BWT.
  @param n_bits Number of bits of each character.
  @param buffer Reduced bwt output. Allocated inside.
*/
int64_t reduce_bwt(char** bwt, uint64_t n, uint n_bits, uint8_t*** buffer, uint steps);

/**
  Function for sorting characters
*/
int char_cmp_func( const void *a, const void *b);

#endif
