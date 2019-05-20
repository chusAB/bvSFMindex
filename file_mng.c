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

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef KNL
#include <hbwmalloc.h>
#endif

#include "mem.h"
#include "file_mng.h"

int
file_to_char(const char * file, char ** data)
{
  FILE *f = fopen(file, "rb");
  if (f == NULL)
  {
    fprintf(stderr, "Cannot open file %s \n", file);
    return -1;
  }

  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  fseek(f, 0, SEEK_SET);

  *data = malloc(fsize + 1);
  uint r = fread(*data, sizeof(char), fsize/sizeof(char), f);
  if (r != (fsize/sizeof(char)))
  {
    fprintf(stderr, "Cannot read data from file %s \n", file);
    return -2;
  }
  fclose(f);

  if((*data)[fsize-1] != '\n' && (*data)[fsize-1] != '\r')
    (*data)[fsize] = 0;
  else while((*data)[fsize-1] == '\n' || (*data)[fsize-1] == '\r')
  {
    (*data)[fsize-1] = 0;
    fsize--;
  }

  return 0;
}

ssize_t
read_seq_from_fasta(FILE * fp, char** line)
{
  size_t len = 0;
  ssize_t read = 0;

  // read header
  read = getline(line, &len, fp);
  if (read == -1) return -1;
  if ((*line)[0] != '>') return -1;

  if ((read = getline(line, &len, fp)) == -1)
  {
     fprintf(stderr, "Error parsing FASTA sequence header\n");
     return -1;
  }

  if ( ((*line)[read-1] != '\n') && ((*line)[read-1] != '\r') )
  {
    (*line)[read] = 0;    // if (len == read) ????
  }
  else while ( ((*line)[read-1] == '\n') || ((*line)[read-1] == '\r') )
  {
    (*line)[read-1] = 0;
    read--;
  }

  return read;
}
