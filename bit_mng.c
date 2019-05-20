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

#include "bit_mng.h"

/* write least significant bits */
void
write_char_to_buffer(uint8_t* buffer, uint n_bits, uint64_t position, uint8_t value)
{
  uint byte_pos = position / 8;
  uint8_t pos_in_byte = position % 8;
  uint8_t value_aux, mask;

  assert(n_bits <= 8);
  /* cross byte border? */
  if (pos_in_byte + n_bits <= 8)
  {
      /* reset bits in the source byte */
      mask = ((1 << n_bits) - 1) << (8 - pos_in_byte - n_bits);
      buffer[byte_pos] = buffer[byte_pos] & (~mask);
      /* write new bits */
      mask = ((1 << n_bits) - 1);
      value = (value & mask) << (8 - pos_in_byte - n_bits);
      buffer[byte_pos] = buffer[byte_pos] | value;
  }
  else
  {
      /* first up to the byte border (8 - pos_in_byte) */
      uint nbits_aux = 8 - pos_in_byte;
      /* reset bits in the source byte */
      mask = (1 << nbits_aux) - 1;
      buffer[byte_pos] = buffer[byte_pos] & (~mask);
      /* write new bits */
      value_aux = value >> (n_bits - nbits_aux);
      buffer[byte_pos] = buffer[byte_pos] | value_aux;
      // printf("mask: 0x%x - value_aux: 0x%x\n", mask, value_aux);

      /* and then from the border of the byte (n_bits - 8 + pos_in_byte) */
      nbits_aux = n_bits - nbits_aux;
      /* reset bits in the source byte */
      mask = ((1 << nbits_aux) - 1) << (8 - nbits_aux);
      buffer[byte_pos+1] = buffer[byte_pos+1] & (~mask);
      /* write new bits */
      mask = ((1 << nbits_aux) - 1);
      value_aux = (value & mask) << (8 - nbits_aux);
      buffer[byte_pos+1] = buffer[byte_pos+1] | value_aux;
      // printf("mask: 0x%x - value_aux: 0x%x\n", mask, value_aux);
  }
}

uint8_t
read_char_from_buffer(uint8_t* buffer, uint n_bits, uint64_t position)
{
  uint64_t byte_pos = position/8;
  uint64_t pos_in_byte = position % 8;
  uint8_t value, value_aux, mask;

  assert(n_bits <= 8);

  /* cross byte border? */
  if (pos_in_byte + n_bits <= 8)
  {
      mask = (1 << n_bits) - 1;
      value = buffer[byte_pos] >> (8 - pos_in_byte - n_bits);
      value = value & mask;
  }
  else
  {
      /* first up to the byte border */
      uint nbits_aux = 8 - pos_in_byte;
      mask = (1 << nbits_aux) - 1;
      value = buffer[byte_pos] & mask;
      // printf("(%x-", value);

      /* and then from the border of the byte */
      nbits_aux = n_bits - nbits_aux;
      mask = (1 << nbits_aux) - 1;
      value_aux = buffer[byte_pos+1] >> (8 - nbits_aux);
      value_aux = value_aux & mask;
      // printf("%x) ", value_aux);

      /* combine the two sequences of bits */
      value = (value << nbits_aux) | value_aux;
  }
  return value;
}
