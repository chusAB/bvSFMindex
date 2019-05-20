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
 
// sched_getcpu()
// #define _GNU_SOURCE

#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <xmmintrin.h>
#include <float.h>
//#include <sched.h>  // sched_getcpu()
#include <sys/syscall.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include "../aux.h"
#include "../mem.h"
#include "../file_mng.h"
#include "../BWT.h"
#include "../bit_mng.h"
#include "../perf.h"
#include "k2d64bv.h"

#ifndef DEBUG_THREADS
    #define DEBUG_THREADS 0
#endif

#ifndef THREADS
    #ifdef KNL
        #define THREADS 256
    #else
        #define THREADS 28
    #endif
#endif

#ifndef NSEQS
    #define NSEQS     4
#endif

#define BYTES_PER_CACHE_BLOCK 64

   
/*  The program runs the kernel nruns times (default 5) and 
 *  reports the *best* result for any iteration after the first,
 *  therefore the minimum value for nruns is 2.
 *  The maximum allowable value for nruns is 128.
 *  Values larger than the default are unlikely to noticeably
 *  increase the reported performance.
 *  nruns can be set on the command line */
#define DEFAULT_RUNS 5

/* disable IACA_START, IACA_END */
#ifndef IACA
    #define IACA_START
    #define IACA_END
#else
    //#include "/opt/intel/iaca-lin64/include/iacaMarks.h"
    #include "/opt/intel/iaca-lin64-v3.0/iacaMarks.h"
#endif

////////////////////////////////////////////////////////////////////////////////
// Module variables
////////////////////////////////////////////////////////////////////////////////

static SFM_t fmi;
static uint64_t mask_64b[64];
static uint32_t nthreads = THREADS;

// input options
static const char *optString = "f:s:t:r:h?";
static const struct option longOpts[] =
{
    {"fmindex",   required_argument,  NULL,   'f'},
    {"sequences", required_argument,  NULL,   's'},
    {"nthreads",  required_argument,  NULL,   't'},
    {"runs",      required_argument,  NULL,   'r'},
    {"help",      no_argument,        NULL,   'h'},
    {NULL,                  0,        NULL,    0 }
};
/*----------------------------------------------------------------------------*/

static struct option_help {
    const char *long_opt, *short_opt, *desc;
} opts_help[] = {
    { "--fmindex", "-f",
      "file storing the fm-index" },
    { "--sequences", "-s",
      "file storing the sequences to search" },
    { "--nthreads", "-t",
      "number of threads" },
    { "--runs", "-r",
      "number of runs" },
    { "--help", "-h",
      "show program usage"},
    { NULL, NULL, NULL }
};

////////////////////////////////////////////////////////////////////////////////
// Private functions
////////////////////////////////////////////////////////////////////////////////

static void
show_usage(char *name, int exit_code)
{
    struct option_help *h;

    printf("usage: %s options\n", name);
    for (h = opts_help; h->long_opt; h++)
    {
        printf(" %s, %s\n ", h->short_opt, h->long_opt);
        printf("    %s\n", h->desc);
    }
    exit(exit_code);
}
/*----------------------------------------------------------------------------*/

// it really computes k2_LF(idx-1)
static __forceinline uint64_t
k2_LF(uint64_t idx, SFM_entry_t *entry)
{
    // uint32_t entry_id     = idx / D_VAL;    // SFM entry index
    uint32_t entry_offset = idx % D_VAL;    // offset en la SFM entry
    uint64_t count = entry->counter;
    count += _popcnt64(entry->data & mask_64b[entry_offset]);
    return count;
}

////////////////////////////////////////////////////////////////////////////////
// Macros
////////////////////////////////////////////////////////////////////////////////

// Check if a sequence has finished the processing
#define _CHECK_FINISHED_SEQ( INDEX )                             \
  if (index[INDEX] < 0 )                                         \
  {                                                              \
    total += end[INDEX] - start[INDEX];                          \
    if (finished_seqs + NSEQS >= block_len)                      \
    {                                                            \
      working_lines[INDEX] = lines[count];                       \
      line_index[INDEX] = count;                                 \
    }                                                            \
    else                                                         \
    {                                                            \
      working_lines[INDEX] = lines_block[finished_seqs + NSEQS]; \
      line_index[INDEX] = bl_offset + finished_seqs + NSEQS;     \
    }                                                            \
    lengths[INDEX] = lines_len[line_index[INDEX]];               \
    finished_seqs++;                                             \
                                                                 \
    uint lut_index = 0, shift_bits = 0;                          \
    uint lut_index_len = 6 - (lengths[INDEX] % 2);               \
    for (uint i=lengths[INDEX] - 1; i > lengths[INDEX] - 1 - lut_index_len; i--)             \
    {                                                                                        \
      lut_index += (uint)(fmi.encoding_table[(uint)working_lines[INDEX][i]]) << shift_bits;  \
      shift_bits += BITS_PER_SYMBOL;                             \
    }                                                            \
    lf += lut_index_len*2;                                             \
    start[INDEX] = fmi.LUT[lengths[INDEX] % KSTEPS][lut_index].start;       \
    end[INDEX]   = fmi.LUT[lengths[INDEX] % KSTEPS][lut_index].end + 1;     \
    index[INDEX] = lengths[INDEX] - KSTEPS - lut_index_len;                 \
    /* printf("\nseq %u: %s\n", line_index[INDEX], working_lines[INDEX]); */ \
    /* decode_symbols(lut_index, seq_tmp, fmi.alphabet, BITS_PER_SYMBOL, lut_index_len); */ \
    /* printf("  LUT_index = %u = %s\n", lut_index, seq_tmp); */             \
  }

// Encode the two next chars to process
#define _ENCODE_CHARS( INDEX )                                             \
  next_symbol[INDEX] = fmi.encoding_table2[                                \
                  *((uint16_t*)(working_lines[INDEX] + index[INDEX]))];    \
  index[INDEX] -= KSTEPS;

static uint64_t __attribute__ ((noinline))
search(char **lines, uint *lines_len, uint count, uint *found, double *glfops)
{
  uint th_bl_size, total = 0;
  uint64_t lf = 0;
  double lfops = 0.0;

  omp_set_num_threads(nthreads);
  th_bl_size = count/nthreads;

  #pragma omp parallel reduction(+:total, lf, lfops) shared(fmi, lines)
  {
    uint64_t start[NSEQS], end[NSEQS];
    SFM_entry_t *start_bl[NSEQS], *end_bl[NSEQS];
    uint8_t next_symbol[NSEQS];
    int index[NSEQS];
    uint line_index[NSEQS], lengths[NSEQS], finished_seqs = 0;
    char * working_lines[NSEQS];
    double start_time, end_time;
    // char seq_tmp[32] = {0};

    // Divide the sequences into blocks, one for each thread
    uint thread_id = omp_get_thread_num();
    uint bl_offset = thread_id*th_bl_size;
    uint block_len  = th_bl_size;
    if (thread_id < count % nthreads)
    {
      block_len++;
      bl_offset += thread_id;
    }
    else
      bl_offset += count % nthreads;

#if DEBUG_THREADS
    uint32_t cpu_num, node_num;
    int status;
    status = syscall(SYS_getcpu, &cpu_num, &node_num, NULL);
    if  (status != -1)
    {
        // int cpu_num = sched_getcpu();
        printf("  th.%2u -> core %2u -> node %u: %9u->%9u\n",
                thread_id, cpu_num, node_num, bl_offset, bl_offset + block_len - 1);
    }
#endif

    char ** lines_block = lines + bl_offset;

    start_time = omp_get_wtime();

    // Get starting positions for the sequences (LUT)
    for(uint j=0; j < NSEQS; j++)
    {
        uint lut_index = 0, lut_index_len = 0, shift_bits = 0;
        working_lines[j] = lines_block[j];
        lengths[j] = lines_len[bl_offset + j];
        line_index[j] = bl_offset + j;
        lut_index_len = 6 - (lengths[j] % 2);
        for (uint i = lengths[j] - 1; i > lengths[j] - 1 - lut_index_len; i--)
        {
            lut_index += (uint)(fmi.encoding_table[(uint)working_lines[j][i]]) << shift_bits;
            shift_bits += BITS_PER_SYMBOL;
        }
        // update stats
        lf += lut_index_len*2;
        start[j] = fmi.LUT[lengths[j] % KSTEPS][lut_index].start;
        end[j]   = fmi.LUT[lengths[j] % KSTEPS][lut_index].end + 1;

        index[j] = lengths[j] - KSTEPS - lut_index_len;
        // Encode *two* chars for the starting pairs (uint16_t)
        next_symbol[j] = fmi.encoding_table2[*((uint16_t*)(working_lines[j] + index[j]))];
        index[j] -= KSTEPS;

#if 0
        printf("\nseq %u: %s\n", line_index[j], working_lines[j]);
        decode_symbols(lut_index, seq_tmp, fmi.alphabet, BITS_PER_SYMBOL, lut_index_len);
        printf("  LUT_index = %u = %s ->", lut_index, seq_tmp);
        printf("  start/end[%d] = %lu/%lu (LUT)\n", j, start[j], end[j]);
        printf("  index[%2u] = %2d/%2d ->", j, index[j] + KSTEPS, lengths[j] - 1);
        decode_symbols(next_symbol[j], seq_tmp, fmi.alphabet, BITS_PER_SYMBOL, KSTEPS);
        printf("  next_symbol[%2d] = %2u = %s ->", j, next_symbol[j], seq_tmp);
        fflush(stdout);
#endif

        // Get starting blocks to search
        start_bl[j] = &fmi.entries[(start[j]/D_VAL)*K2_SYMBOLS + next_symbol[j]];
        end_bl[j]   = &fmi.entries[(  end[j]/D_VAL)*K2_SYMBOLS + next_symbol[j]];

        // Prefetch blocks for next execution of sequence j into L2
        _mm_prefetch((char*) start_bl[j], PREFETCH_HINT_L2);
        _mm_prefetch((char*) end_bl[j],   PREFETCH_HINT_L2);
    }

    while(finished_seqs < block_len)
    {
        /* start of loop to analyze */
        IACA_START

        for (uint j=0; j < NSEQS; j++)
        {
            // Prefetch blocks for sequence j into L1
            #ifdef DUAL_PREFETCH
              _mm_prefetch((char*) start_bl[j], PREFETCH_HINT_L1);
              _mm_prefetch((char*) end_bl[j],   PREFETCH_HINT_L1);
            #endif

            // LFs for sequence j
            start[j] = k2_LF(start[j], start_bl[j]);
            end[j]   = k2_LF(end[j]  , end_bl[j]);

            // printf("  start/end[%2u] = %2lu/%2lu\n", j, start[j], end[j]);
        
            // Check if finished sequence j
            _CHECK_FINISHED_SEQ(j);
            _ENCODE_CHARS(j);

#if 0
            printf("  index[%2u] = %2d/%2d ->", j, index[j] + KSTEPS, lengths[j] - 1);
            decode_symbols(next_symbol[j], seq_tmp, fmi.alphabet, BITS_PER_SYMBOL, KSTEPS);
            printf("  next_symbol[%2d] = %2u = %s ->", j, next_symbol[j], seq_tmp);
            fflush(stdout);
#endif

            // Calculate blocks for sequence j
            start_bl[j] = &fmi.entries[(start[j]/D_VAL)*K2_SYMBOLS + next_symbol[j]];
            end_bl[j]   = &fmi.entries[(end[j]/D_VAL)*K2_SYMBOLS   + next_symbol[j]];

            // Prefetch blocks for next execution of sequence j into L2
            _mm_prefetch((char*) start_bl[j], PREFETCH_HINT_L2);
            _mm_prefetch((char*) end_bl[j],   PREFETCH_HINT_L2);
            ////////////////////////////////////////////////////////////////////
        }
        // update stats
        lf += 2*KSTEPS*NSEQS;
    }

    end_time = omp_get_wtime();
    lfops = lf/(end_time - start_time);

    /* end of loop to analyze */
    IACA_END

    #pragma omp barrier

  } //#pragma omp parallel

  (*found) = total;
  (*glfops) = lfops/GIGA;
  return lf;
}

static void
metrics(double *sample, uint64_t *lf, double *sample_glfops, int nruns)
{
    // Compute results
    double time_min = FLT_MAX, time_max = 0;
    double glfops_min = FLT_MAX, glfops_max = 0;
    double time_sum = 0, time_deviation_sum = 0, time_avg, time_deviation;
    uint32_t time_min_idx, time_max_idx;
    double glfops_sum = 0.0, glfops_avg;
    double GLF = lf[0]/GIGA;

    if (!elements_64b_equal(lf, nruns))
    {
        printf("ERROR: different runs performed different number of LFOPS\n");
        return;
    }

    for (int i = 1; i < nruns; i++)
    {
        time_sum += sample[i];
        if (sample[i] < time_min)
        {
            time_min = sample[i];
            time_min_idx = i;
        }
        if (sample[i] > time_max)
        {
            time_max = sample[i];
            time_max_idx = i;
        }

        glfops_sum += (sample_glfops[i]);
        if (sample_glfops[i] < glfops_min)
        {
            glfops_min = sample_glfops[i];
        }
        if (sample_glfops[i] > glfops_max)
        {
            glfops_max = sample_glfops[i];
        }
    }
    
    time_avg = time_sum / (nruns - 1);
    glfops_avg = glfops_sum / (nruns - 1);
    
    for (int i = 1; i < nruns; i++)
    {
        time_deviation_sum += (sample[i] - time_avg) * (sample[i] - time_avg);
    }
    time_deviation = sqrt(time_deviation_sum / (nruns-1));

    // Show metrics
    printf("Best throughput: %6.3f GLFOPS (%.3f)\n", GLF/time_min, glfops_max);
    printf("Avg. throughput: %6.3f GLFOPS (%.3f)\n", GLF/time_avg, glfops_avg);
    printf("Throughputs: [%.3f]", GLF/sample[0]);
    for (int i = 1; i < nruns; i++)
    {
        printf(" %.3f", GLF/sample[i]);
    }
    printf("\n");    
    printf("Min time: %6.3f s (run #%u)\n", time_min, time_min_idx);
    printf("Avg time: %6.3f s\n", time_avg);
    printf("Max time: %6.3f s (run #%u)\n", time_max, time_max_idx);    
    printf("Times: [%.3f]", sample[0]);
    for (int i = 1; i < nruns; i++)
    {
        printf(" %.3f", sample[i]);
    }
    printf("\n");    
    printf("Standard deviation = %.3f%%\n", (time_deviation / time_avg) * 100);
}

////////////////////////////////////////////////////////////////////////////////
// Main function
////////////////////////////////////////////////////////////////////////////////

int
main(int argc, char *argv[])
{
  uint total[MAXRUNS], count = 0, lines_size = 1000;
  int nruns = DEFAULT_RUNS;
  FILE * fp;
  char ** lines;
  uint * lines_len;
  ssize_t read;
  double start_timer, end_timer, sample[MAXRUNS];
  double start0, end0;
  double sample_glfops[MAXRUNS];
  uint64_t bases = 0, lf[MAXRUNS] = { 0 };
  char *fmi_file = 0;
  char *seq_file = 0;
  int n = 0, option = 0;

  printf("Program version: 20190206\n");
  printf(HLINE);

  while(1)
  {
      option = getopt_long(argc, argv, optString, longOpts, NULL /* &longIndex */);
      if (option == -1) break;
    
      switch(option)
      {
          case 'f':
              fmi_file = optarg;
              break;

          case 's':
              seq_file = optarg;
              break;
    
          case 't':
              n = sscanf(optarg, "%u", &nthreads);
              if ((n == -1) || (nthreads < 1))
              {
                  printf("ERROR: wrong number of threads\n\n");
                  exit(1);
              }
              break;

          case 'r':
              n = sscanf(optarg, "%u", &nruns);
              if ((nruns < 2) || (nruns > MAXRUNS))
              {
                  printf(HLINE);
                  printf("** Unsupported number of runs: changed to %d (default) **\n", DEFAULT_RUNS);
                  nruns = DEFAULT_RUNS;
              }
              break;

          case 'h':
              show_usage(argv[0], 0);
              break;
    
          default:
              show_usage(argv[0], 1);
      }
  }

  /* check arguments */
  if (fmi_file == 0)
  {
      printf("ERROR: fm-index file not specified\n");
      show_usage(argv[0], 1);
  }
  if (seq_file == 0)
  {
      printf("ERROR: sequences file not specified\n");
      show_usage(argv[0], 1);
  }

#if LIBNUMA
#ifndef KNL
  // memory policy configuration
  mem_conf();
  printf(HLINE);
#endif
#endif

  printf("Loading FM-index...\n");
  start_timer = omp_get_wtime();
  load_SFM(fmi_file, &fmi);
  end_timer = omp_get_wtime();
  printf("OK. Index loaded in %fs\n", end_timer - start_timer);

#if LIBNUMA
  // free huge pages in each node
  hugepages_status();
  printf(HLINE);
#endif

  // Loading sequences into memory
  printf("Loading sequences...\n");
  start_timer = omp_get_wtime();
  lines = malloc(lines_size*sizeof(char*));
  lines_len = malloc(lines_size*sizeof(uint));
  if ( (lines == NULL) || (lines_len == NULL) )
  {
    printf("Error at malloc\n");
    exit(EXIT_FAILURE);
  }

  fp = fopen(seq_file, "r");
  if (fp == NULL)
  {
    printf("Error opening file %s\n", seq_file);
    exit(EXIT_FAILURE);
  }

  while ((read = read_seq_from_fasta(fp, &(lines[count]))) >= 0)
  {
    lines_len[count] = read;
    bases += read;
    count++;
    if (count >= lines_size)
    {
      lines_size += 1000;
      lines = realloc(lines, lines_size*sizeof(char*));
      lines_len = realloc(lines_len, lines_size*sizeof(uint));
    }
  }
  fclose(fp);
  end_timer = omp_get_wtime();
  printf("OK. %.2f Msequences loaded in %fs (%.3f Mseq/s)\n",
          count/MEGA, end_timer - start_timer, (double)(count)/(MEGA*(end_timer - start_timer)));
  printf(HLINE);
  fflush(stdout);

  lines = realloc(lines, (count+1)*sizeof(char*));
  lines_len = realloc(lines_len, (count+1)*sizeof(uint));
  lines[count] = fmi.start;
  lines_len[count] = strlen(fmi.start);

  printf("Parameters\n");
  printf("- FM-index file: %s\n", fmi_file);
  printf("- Index size: %.1fGiB (%lu characters)\n", (double)(fmi.len)/GiB, fmi.len);
  printf("- Number of threads: %d\n", nthreads);
  printf("- Overlapped sequences: %d\n", NSEQS);
  printf("- Sequence file: %s\n", seq_file);
  printf("- Number of bases: %.2f Gbases (%lu)\n", bases/GIGA, bases);
  printf("- Number of sequences: %.2f Mseq (%u)\n", count/MEGA, count);
  printf("- Average sequence length: %.2f bases\n", (double) bases/count);
  printf("- Number of sequences processed per thread: %.1fM (%u)\n", (double) (count)/(MEGA*nthreads), count/nthreads);
  printf(HLINE);
  printf("The kernel will be executed %d times.\n", nruns);
  printf("The *best* throughput will be reported (excluding the first iteration).\n");
  printf(HLINE);

  // Mask initialization
  mask_init(mask_64b);

  // Executing the FM-index count
  printf("Starting search... \n");
  
#ifdef PERF
  perf_init();
#endif

  start0 = omp_get_wtime();

  for (int run = 0; run < nruns; run++)
  {
#ifdef PERF
      perf_enable();
#endif

      start_timer = omp_get_wtime();
      lf[run] = search(lines, lines_len, count, &total[run], &sample_glfops[run]);
      end_timer = omp_get_wtime();
      sample[run] = end_timer - start_timer;

#ifdef PERF
      perf_stop();
      perf_read_sample(run);
      perf_reset();
#endif
  }

  end0 = omp_get_wtime();
  printf("OK\n");

  /* Expected bases do not match with processed bases because fmi.start characters are processed */
  /* if (lf != bases*2) printf("Error counting LF(): expected %lu but counted %lu\n", bases*2, lf); */

  printf("Occurrences found:");
  if (elements_equal(total, nruns))
  {
      printf(" %u (in each run)\n", total[0]);
  }
  else
  {
      for (int i = 0; i < nruns; i++)
          printf(" %u", total[i]);
      printf("\n");
  }

#if 0
  printf("Total processed bases: %.2fG (expected %.2fG = nbases x nruns x 2 lfs/base)\n", sum(lf, nruns)/(2.0*GIGA), nruns*bases/GIGA);
  printf("  ");
  for (int i = 0; i < nruns; i++)
      printf("%.2fG (%lu) ", lf[i]/(2.0*GIGA), lf[i]/2);
  printf("\n");
#endif
  printf("Total LFOP: %.2fG (expected %.2fG)\n", sum(lf, nruns)/GIGA, 2*nruns*bases/GIGA);
  printf("Total time: %f\n", end0 - start0);
  printf("Raw throughput: %6.3f GLFOPS\n", sum(lf, nruns)/(end0 - start0)/GIGA);
  printf(HLINE);

  metrics(sample, lf, sample_glfops, nruns);
  printf(HLINE);

#ifdef PERF
  perf_print_samples(nruns);
  perf_close();
  printf(HLINE);
#endif

  if (lines)
  {
    for(uint i=0; i<count; i++) free(lines[i]);
    free(lines);
  }
  free_SFM(&fmi);
  return 0;
}
