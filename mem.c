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
#include <assert.h>
#include <string.h>
#include <sched.h> // For sched_setaffinity
#include <unistd.h>
#if LIBNUMA
#include <numa.h>
#include <numaif.h>

#include "mem.h"


static void
print_oneline_file(char *filename)
{
  size_t len = 0;
  FILE * fp;
  char *linep = NULL;
  
  fp = fopen(filename, "r");
  if (fp == NULL)
  {
    printf("    NA");
    return;
  }

  if (-1 == getline(&linep, &len, fp))
	  printf("    NA");
  else
  {
      size_t linelen = strlen(linep);
      while ( (linep[linelen-1] == '\n') || (linep[linelen-1] == '\r') )
      {
        linep[linelen-1] = 0;
        linelen--;
      }
      printf("%6s", linep);
  }

  free(linep);
  fclose(fp);
}


void
hugepages_status()
{
  // Get numa configuration
  int numa_nodes;
  char string[512];

  // check if system supports NUMA API
  int available = numa_available();
  if (available == -1)
  {
      numa_nodes = 0;
      return;
  }
  else
      numa_nodes = numa_num_configured_nodes();

  printf("Free hugepages   2MB   1GB\n");
  for (int j=0; j < numa_nodes; j++)
  {
    printf("node %u:       ", j);
    sprintf(string, "/sys/devices/system/node/node%u/hugepages/hugepages-2048kB/free_hugepages", j);
    print_oneline_file(string);

    sprintf(string, "/sys/devices/system/node/node%u/hugepages/hugepages-1048576kB/free_hugepages", j);
    print_oneline_file(string);

    printf("\n");
  }
  // printf("\n\n");
}


int
mem_conf()
{
  // uint64_t nodes_mask = 0;
  int ncores = sysconf(_SC_NPROCESSORS_ONLN);
  int core = sched_getcpu();
  int node = 0, max_node = 0;
  int available = numa_available();

  // check if system supports NUMA API
  if (available == -1) return 0;

  max_node = numa_max_node();
  printf("System info: %d %s, %d cores (%d cores/node)\n",
          max_node+1, max_node == 0? "node":"nodes", ncores, ncores/(max_node+1));

  node = numa_node_of_cpu(core);
  printf("Process info: running in node %d, core %d.\n", node, core);

  return node;
}
#endif
