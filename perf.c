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

#define _GNU_SOURCE 1

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
#include <time.h>
#include <errno.h>
#include "perf.h"

/* static variables */

static char *events[HW_CTRS] = {"Instructions", "Branches", "Branch Misses", "DTLB Read Misses"};
static char *short_events[HW_CTRS] = {"Ir", "Bc", "Bcm", "DTLBm"};  // valgrind notation
static int fd[HW_CTRS];
static uint64_t hw_counters[MAXRUNS][HW_CTRS] = {{0}};

/* forward declarations */
void perf_reset();


/* ************************************************************************** */
/* static functions */
/* ************************************************************************** */

// Open a perf event
static int
perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                int cpu, int group_fd, unsigned long flags)
{
    return syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
}

// Generate the perf struct
static void
perfStruct(struct perf_event_attr *pe, int type, int config)
{
    memset(pe, 0, sizeof(struct perf_event_attr));
    pe->size = sizeof(struct perf_event_attr);
    pe->disabled = 1;

    // Exclude kernel and hipervisor from being measured
    pe->exclude_kernel = 1;
    pe->exclude_hv = 1;

    /* children inherit it */
    pe->inherit=1;

    // Type of event to measure
    pe->type = type;
    pe->config = config;
}

/* eases parsing */
static void
perf_print_one_line(uint64_t *count)
{
    // print short events
    printf("\n");
    for (int i = 0; i < HW_CTRS; i++) printf("  %14s", short_events[i]);
    printf("  %8s\n", "Bcmpki");

    // print counters
    for (int i = 0; i < HW_CTRS; i++) printf("  %14lu", count[i]);

    /* Branch conditional prediction misses per kiloinstruction */
    printf("  %8.2f", (1000.0*count[2])/count[0]);
    printf("    hw_counters\n");
}

/* print counters values stored in the variable passed as argument */
static void
perf_print_counters(uint64_t *count, char *title)
{
    printf("%s\n", title);
    for (int i = 0; i < HW_CTRS; i++)
        printf("  %16s: %10.2fM\n", events[i], count[i]/1e6);

    /* Branch conditional prediction misses per kiloinstruction */
    printf("  %16s: %10.2f\n", "Bcmpki", (1000.0*count[2])/count[0]);

    /* eases parsing */
    perf_print_one_line(count);
}


/* ************************************************************************** */
/* public functions */
/* ************************************************************************** */

/* initializes hardware counters */
int
perf_init()
{
    struct perf_event_attr pe[HW_CTRS];
    // Type of events
    int type[HW_CTRS] = {PERF_TYPE_HARDWARE, PERF_TYPE_HARDWARE, PERF_TYPE_HARDWARE, PERF_TYPE_HW_CACHE};
    // Name of the events to measure
    int event[HW_CTRS] = {PERF_COUNT_HW_INSTRUCTIONS, PERF_COUNT_HW_BRANCH_INSTRUCTIONS, PERF_COUNT_HW_BRANCH_MISSES, PERF_COUNT_HW_CACHE_DTLB | (PERF_COUNT_HW_CACHE_OP_READ << 8) | (PERF_COUNT_HW_CACHE_RESULT_MISS << 16)};

    // Get perfStruct
    for (int i = 0; i < HW_CTRS; ++i) perfStruct(&pe[i], type[i], event[i]);
    // Open perf leader
    for (int i = 0; i < HW_CTRS; ++i)
    {
        if ((fd[i] = perf_event_open(&pe[i], 0, -1, -1, 0)) == -1)
        {
            fprintf(stderr,"Error opening event %llx %s\n", pe[i].config, strerror(errno));
            return -1;
        }
    }
    perf_reset();
    return 0;
}

/* resets hardware counters */
void
perf_reset()
{
    // Enable descriptor to read hw counters
    for (int i = 0; i < HW_CTRS; ++i)
    {
        ioctl(fd[i], PERF_EVENT_IOC_RESET, 0);
    }
}

/* activates hardware counters */
void
perf_enable()
{
    // Enable descriptor to read hw counters
    for (int i = 0; i < HW_CTRS; ++i)
    {
        ioctl(fd[i], PERF_EVENT_IOC_ENABLE, 0);
    }
}

/* initializes, resets and activates hardware counters */
int
perf_start()
{
    int n;

    n = perf_init();
    if (n == -1) return -1;
    perf_reset();
    perf_enable();
    return 0;
}

/* deactivates hardware counters (stops counting events) */
void
perf_stop()
{
    // "Turn off" hw counters
    for (int i = 0; i < HW_CTRS; ++i) ioctl(fd[i], PERF_EVENT_IOC_DISABLE, 0);
}

/* read hardware counters */
int
perf_read_sample(uint32_t run)
{
    // Reading HW counters
    for (int i = 0; i < HW_CTRS; i++)
    {
        if (sizeof(uint64_t) != read(fd[i], &hw_counters[run][i], sizeof(uint64_t)))
        {
            printf("\nERROR: lectura incorrecta en contador %16s\n", events[i]);
            return -1;
        }
    }
    return 0;
}

/* closes the descriptors associated with the counters */
void
perf_close()
{
    // Close file descriptors
    for (int i = 0; i < HW_CTRS; ++i) close(fd[i]);
}

/* prints counter values and closes descriptors */
void
perf_print()
{
    // Reading HW counters
    perf_read_sample(0);

    // Close file descriptors
    perf_close();

    // print counters
    perf_print_counters(hw_counters[0], "Hardware counters");
}

void
perf_print_samples(int nruns)
{
    uint64_t count_avg[HW_CTRS] = { 0 };
 
    printf("Hardware counters samples\n");
    printf("      \t");
    for (int j = 0; j < HW_CTRS; j++)
        printf("%14s\t", short_events[j]);
    printf("\n");

    // Print samples and compute average
    for (int i = 0; i < nruns; i++)
    {
        printf("run #%d\t", i);
        for (int j = 0; j < HW_CTRS; j++)
        {
            if (i != 0) count_avg[j] += hw_counters[i][j];
            printf("%14lu\t", hw_counters[i][j]);
        }
        printf("\n");
    }
    printf("avg.  \t");
    for (int j = 0; j < HW_CTRS; j++)
    {
        count_avg[j] /= (nruns-1);
        printf("%14lu\t", count_avg[j]);
    }
    printf("\n\n");

    perf_print_counters(count_avg, "Hardware counters average");
}
