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

#include <assert.h>
#include <stdio.h>

#include "aux.h"

/* Fast ceiling of an unsigned integer division */
uint64_t
ceil_uint_div(uint64_t x, uint64_t y)
{
    return (x/y + (x % y != 0));
}


int
elements_64b_equal(uint64_t *v, int vlen)
{
  for (int i=1; i < vlen; i++)
  {
    if (v[i] != v[0])
    {
      printf("run #%d: %lu vs. run #0: %lu\n", i, v[i], v[0]);
      return 0;
    }
  }
  return 1;
}

int
elements_equal(uint32_t *v, int vlen)
{
  for (int i=1; i < vlen; i++)
  {
    if (v[i] != v[0])
    {
      printf("run #%d: %u vs. run #0: %u\n", i, v[i], v[0]);
      return 0;
    }
  }
  return 1;
}

uint64_t
sum(uint64_t *v, int vlen)
{
  uint64_t vsum = 0;

  for (int i=0; i < vlen; i++)
     vsum += v[i];

  return vsum;
}
