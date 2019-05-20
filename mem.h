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

#ifndef _MEM_H_
#define _MEM_

#include <stdio.h>
#include "types.h"

/* For explicit huge page allocation */
#define ADDR             (void *)(0x0UL)
#define PROTECTION       (PROT_READ | PROT_WRITE)

#define HUGE_PAGE_FLAGS  (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB)
// #define MAP_HUGE_SHIFT   26
// #define MAP_HUGE_1GB     (30 << MAP_HUGE_SHIFT)
// #define HUGE_PAGE_FLAGS  (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_HUGE_1GB)

#define ALIGN_2MB 1U << 21

#define V_1MB    (1024.0*1024.0)
#define V_1GB    (V_1MB*1024.0)

#define PREFETCH_HINT_1ST _MM_HINT_T1
#define PREFETCH_HINT_2ND _MM_HINT_T0
#define PREFETCH_HINT_L2 _MM_HINT_T1
#define PREFETCH_HINT_L1 _MM_HINT_T0

void hugepages_status();

int mem_conf();

#endif
