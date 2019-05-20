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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <divsufsort64.h>

#include "BWT.h"
#include "aux.h"
#include "bit_mng.h"
// #include "divsufsort64.h"

int
get_unique_elements(char * in, char ** out, uint n)
{
  uint i;
  int j=0, size=20, count=0;

  *out = calloc(20, sizeof(char));
  // *out = (char*)malloc(20*sizeof(char));
  if (*out == NULL)
  {
    fprintf(stderr, "Error at output malloc for get_unique_elements \n");
    return -1;
  }

  for(i = 0; i < n; i++)
  {
    if (in[i] != '$')
    {
      for (j=0; j<count; j++)
      {
        if (in[i] == (*out)[j])
         break;
      }
      if (j == count)
      {
        if (count >= size)
        {
          *out = (char*) realloc(*out, (size+20)*sizeof(char));
          size += 20;
        }
        (*out)[count] = in[i];
        count++;
      }
    }
  }
  // size adjustment
  *out = (char*)realloc(*out, count*sizeof(char));
  return count;
}

/* in: input string
 * out: BWT
 * n: length of input string */
int
get_bwt(char** in, char*** out, uint64_t n, uint steps, uint64_t** end)
{
  saidx64_t * SA;

  // We add $ at the end of the string
  *in = realloc(*in, (n+2)*sizeof(char));
  (*in)[n] = '$'; (*in)[n+1] = 0;  //strcat(*in, "$");
  // before: n chars, n+1 allocated bytes
  // after: n+1 chars, n+2 allocated bytes

  // Memory allocation for the suffix array (SA)
  SA = (saidx64_t*) malloc((n+1)*sizeof(saidx64_t));
  if (SA == NULL)
  {
    fprintf(stderr, "Error at malloc for SA: %lu bytes (%.2f GB) requested but not allocated\n", (n+1)*sizeof(saidx64_t), (float) (n+1)*sizeof(saidx64_t)/1.0e9);
    return -1;
  }

  *end = (uint64_t*) malloc(steps*sizeof(uint64_t));
  if (*end == NULL)
  {
    fprintf(stderr, "Error at malloc for $ pos \n");
    return -2;
  }

  // suffix array calculation
  if (divsufsort64((sauchar_t*)(*in), SA, n+1) < 0)
  {
    fprintf(stderr, "Error when generating SA \n");
    return -3;
  }

  // Out memory allocation
  *out = (char**) malloc(steps*sizeof(char*));
  if (*out == NULL)
  {
    fprintf(stderr, "Error at malloc for BWT output \n");
    return -4;
  }
  for(uint64_t i = 0; i < steps; i++)
  {
    (*out)[i] = (char*) malloc(sizeof(char)*(n+1));
    if ((*out)[i] == NULL)
    {
      fprintf(stderr, "Error at malloc for BWT output \n");
      return -4;
    }
  }

  // BWT calculation from SA
  // for (uint64_t i=0; i<n+1; i++) {
  //    if (SA[i] > 0)      BWT[i] = T[SA[i]-1];
  //    else                BWT[i] = '$'; }
  for(uint64_t i=0; i < n+1; i++)
  {
    for(int64_t j=0; j < steps; j++)
    {
      if (SA[i] > j)
      {
        (*out)[j][i] = (*in)[SA[i]-1-j];
      }
      else if (SA[i] < j)
      {
        (*out)[j][i] = (*in)[n-j+SA[i]];
      }
      else
      {
        (*end)[j] = i;
        (*out)[j][i] = '$';  // (*out)[j][i] = (*in)[n];
      }
    }
  }

  // 0 ended string
  for (int64_t i=0; i < steps; i++)
    (*out)[i][n+1] = 0;

  free(SA);
  return 0;
}

int
encode_bwt(char** bwt, char* codes, uint n_bwt, uint n_chars, uint steps, uint64_t* end)
{
  for(uint64_t k = 0; k < steps; k++)
  {
    for(uint64_t i = 0; i < n_bwt; i++)
    {
      for(uint64_t j = 0; j < n_chars; j++)
      {
        if (bwt[k][i] == codes[j])
        {
          bwt[k][i] = j;
          break;
        }
      }
    }
    // $ coding
    bwt[k][end[k]] = 0;
  }
  return 0;
}


int64_t
reduce_bwt(char** bwt, uint64_t n, uint n_bits, uint8_t*** buffer, uint steps)
{
  uint64_t i,j;
  uint64_t n_bytes = ceil_uint_div(n*(uint64_t)n_bits, 8);

  *buffer = (uint8_t**) malloc(steps*sizeof(uint8_t*));
  if (*buffer == NULL)
  {
    fprintf(stderr, "Error at malloc for BWT output \n");
    return -4;
  }
  for(i = 0; i<steps; i++)
  {
    (*buffer)[i] = calloc(n_bytes, sizeof(uint8_t));
    if ((*buffer)[i] == NULL)
    {
      fprintf(stderr, "Error at calloc for BWT buffer \n");
      return -4;
    }
  }

  for(j=0; j<steps; j++)
  {
    for(i=0; i<n; i++)
      write_char_to_buffer((*buffer)[j], n_bits, ((uint64_t)n_bits)*i, bwt[j][i]);
  }

  return n_bytes;
}


int
char_cmp_func(const void *a, const void *b)
{
  return *(char*)a - *(char*)b;
}
