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
 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#include "../file_mng.h"
#include "../BWT.h"
#include "../bit_mng.h"
#include "k2d64bv.h"

/* return wall time in seconds */
static double
get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time,NULL)) {
        exit(-1); // return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int
main(int argc, const char *argv[])
{
  char * data;
  uint64_t data_len;
  char ** bwt;
  char * unique_data;
  char outfile[80];
  uint8_t ** reduced;
  SFM_t fmi;
  uint n_bits;
  uint64_t *end;
  int64_t reduced_len;
  uint64_t C[KSTEPS][SYMBOLS];
  double wall_0, wall_1;
  int n, unique_len;
  int verbose = 0;

  if (argc < 2)
  {
    fprintf(stderr, "Use: ./%s file\n", argv[0]);
    return 0;
  }
  if (argc == 3)
    verbose = 1;

  printf("Reading FM-index file %s... ", argv[1]);
  wall_0 = get_wall_time();
  n = file_to_char(argv[1], &data);
  if (n < 0) exit(1);
  data_len = strlen(data);
  wall_1 = get_wall_time();
  printf("OK\n");
  printf("Total time: %.3fs\n", wall_1 - wall_0);
  if (verbose)
  {
      printf("Reference text");
      dump_array(data, data_len);
  }
  /*--------------------------------------------------------------------------*/

  // init_C(C);
  // build_C(&C[0][0], data, data_len);
  // printf("---C Array---(%lu elements)\n", data_len);
  // dump_C(C);
  /*--------------------------------------------------------------------------*/

  printf("Getting BWT of %lu characters... \n", data_len);
  wall_0 = get_wall_time();
  n = get_bwt(&data, &bwt, data_len, KSTEPS, &end);
  if (n < 0) exit(1);
  data_len++;  // $ character
  wall_1 = get_wall_time();
  printf("OK\nBWT Generated. Length: %lu\n", data_len);
  printf("Total time: %.3fs\n", wall_1 - wall_0);
  if (verbose)
  {
      dump_BWT(bwt, KSTEPS, data_len);
      for (int i = 0; i < KSTEPS; i++) printf("- end[%d] = %lu\n", i, end[i]);
  }
  /*--------------------------------------------------------------------------*/

  // init_C(C);
  // for(int j=0; j<KSTEPS; j++)
  //  build_C(&C[j][0], &bwt[j][0], data_len);

  // printf("---C Array---(%lu elements)\n", data_len);
  // dump_C(C);
  /*--------------------------------------------------------------------------*/

  printf("Encoding BWT...\n");
  wall_0 = get_wall_time();
  unique_len = get_unique_elements(bwt[0], &unique_data, data_len);
  if (unique_len < 0) exit(1);

  n_bits = (uint)(ceil(log2(unique_len)));
  printf(" -> %d Unique elements. Each character can be stored in %u bits\n", unique_len, n_bits);
  // dump_array(unique_data, unique_len);
  qsort((void*)unique_data, unique_len, sizeof(char), char_cmp_func);
  dump_array(unique_data, unique_len);
  // printf("unique_data: %s\n", unique_data);

  encode_bwt(bwt, unique_data, data_len, unique_len, KSTEPS, end);
  wall_1 = get_wall_time();
  printf("OK, encoded BWT\n");
  printf("Total time: %.3fs\n", wall_1 - wall_0);
  if (verbose)
  {
      dump_encoded_BWT(bwt, KSTEPS, data_len, unique_data, end);
      for (int i = 0; i < KSTEPS; i++) printf("- end[%d] = %lu\n", i, end[i]);
  }
  /*--------------------------------------------------------------------------*/

  init_C(C);
  for(int j=0; j < KSTEPS; j++)
  {
    for(uint64_t i=0; i<data_len; i++)
      C[j][(uint)bwt[j][i]]++;
  }
  C[0][0]--;  // $ se codifica como 0
  C[1][0]--;  // $ se codifica como 0
  printf("---Encoded C Array---(%lu elements)\n", data_len);
  dump_C(C);
  /*--------------------------------------------------------------------------*/

  printf("Compressing BWT... ");
  wall_0 = get_wall_time();
  reduced_len = reduce_bwt(bwt, data_len, n_bits, &reduced, KSTEPS);
  if (reduced_len < 0) exit(1);
  wall_1 = get_wall_time();
  printf("OK\n -> Compression ratio: %.2f = %lu / %lu (original/compressed bytes)\n",
          (float) data_len / (reduced_len*KSTEPS), data_len, reduced_len*KSTEPS);
  printf("Total time: %.3fs\n", wall_1 - wall_0);
  free(bwt);
  /*--------------------------------------------------------------------------*/

  printf("Generating FM-index... ");
  wall_0 = get_wall_time();
  fmi.start = malloc(sizeof(char)*501);
  memcpy(fmi.start, data, 500);
  fmi.start[500] = 0;
  free(data);
  // dump_array(unique_data, unique_len);
  generate_SFM(&fmi, reduced, data_len, unique_data, 64, end);

  if (verbose) dump_SFM(&fmi);

  snprintf(outfile, sizeof(outfile), "%s.k%dd%dbv.fmi", argv[1], KSTEPS, D_VAL);
  write_SFM(outfile, &fmi);
  wall_1 = get_wall_time();
  printf("OK -> FM-index written to file %s\n", outfile);
  printf("FM-index time: %.3fs\n", wall_1 - wall_0);
  printf("-------------------------------------------------\n\n");

  /*--------------------------------------------------------------------------*/

  exit(0);
}
