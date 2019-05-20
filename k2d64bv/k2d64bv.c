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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <ctype.h>
#include <sys/mman.h>
#ifdef KNL
#include <hbwmalloc.h>
#endif
//#include <divsufsort.h>

#include "../aux.h"
#include "../mem.h"
#include "../file_mng.h"
#include "../BWT.h"
#include "../bit_mng.h"
#include "k2d64bv.h"


////////////////////////////////////////////////////////////////////////////////
// Module variables

static SFM_t *pfmi;
static SFM_entry_t *ROcc;
static uint64_t mask_64b[64];

static __forceinline uint64_t
k2_LF(unsigned char symbol, uint64_t idx)
{
    uint32_t entry_id        = idx / D_VAL;     // SFM entry index
    uint32_t entry_offset    = idx % D_VAL;     // offset en la SFM entry
    SFM_entry_t entry = ROcc[entry_id*K2_SYMBOLS + symbol];
    uint64_t count = entry.counter;

    //printf("Count A: %lu %lu %lu\n", count, entry.data, idx);
    count += _popcnt64(entry.data & mask_64b[entry_offset]);
    //printf("Count B: %lu %lu \n", count, mask_64b[entry_offset]);
    return count;
}

void
init_C(uint64_t C[KSTEPS][SYMBOLS])
{
  for(int i = 0; i < KSTEPS; i++)
  {
    for(int j = 0; j < SYMBOLS; j++)
      C[i][j] = 0;
  }
}

void
dump_C(uint64_t C[KSTEPS][SYMBOLS])
{
  for(int i=0; i < KSTEPS; i++)
  {
    printf("C[%i]: ", i);
    for(int j=0; j < SYMBOLS; j++)
      printf(" %lu ", C[i][j]);
    printf("\n");
  }
}

void
dump_C_SFM(uint64_t C[K2_SYMBOLS+1])
{
  for(int i=0; i < K2_SYMBOLS + 1; i++)
    printf("C[%i]: %lu ", i, C[i]);

  printf("\n");
}

/* nrows: KSTEPS */
/* matrix[KSTEPS][LEN] */
void
dump_BWT(char **matrix, uint nrows, uint ncols)
{
  printf("BWT\n  ");
  for(uint i=0; i < D_VAL; i++)
  {
    printf("%2u ", i);
  }
  printf("\n  ");

  for(uint i = 0; i < ncols; i++)
  {
    for(int j = nrows - 1; j >= 0; j--)
    {
      printf("%c", matrix[j][i]);
    }
    printf(" ");
    if (i % D_VAL == D_VAL - 1) printf("\n  ");
  }
  printf("\n");
}

/* nrows: KSTEPS */
/* matrix[KSTEPS][LEN] */
void
dump_encoded_BWT(char **matrix, uint nrows, uint ncols, char *codes, uint64_t* end)
{
  printf("BWT\n  ");
  for(uint i=0; i < ncols; i++)
  {
    for(int j = nrows - 1; j >= 0; j--)
    {
      if (i == end[j]) printf("$");
      else             printf("%c", codes[(unsigned char) matrix[j][i]]);
    }
    printf(" ");
    if (i % D_VAL == D_VAL - 1) printf("\n  ");
  }
  printf("\n");
}

#if 0
/* nrows: KSTEPS */
void
dump_encoded_BWT(char **matrix, uint nrows, uint ncols, char *codes, uint64_t* end)
{
  for(uint i=0; i < nrows; i++)
  {
    printf("BWT[%u]\n  ", i);
    for(uint j=0; j < ncols; j++)
    {
      if (j == end[i]) printf("$");
      else             printf("%c", codes[matrix[i][j]]);
      if (j % D_VAL == D_VAL - 1) printf("\n  ");
    }
    printf("\n");
  }
}
#endif

void
build_C(uint64_t *fc, char *data, uint64_t data_len)
{
  char c;

  for(uint64_t i=0; i < data_len; i++)
  {
    c = tolower(data[i]);
    if (c == 'a')       fc[0]++;
    else if (c == 'c')  fc[1]++;
    else if (c == 'g')  fc[2]++;
    else if (c == 't')  fc[3]++;
  }
}

void
dump_array(char *array, uint64_t len)
{
  printf(" -> ");
  for(uint64_t i=0; i<len; i++)
    printf("%c", array[i]);   // printf("%x|", array[i]);
  printf("\n");
}

void
mask_init(uint64_t *mask)
{
    // it is only necessary to initialize even indexes
    // because offset measured in bits will always be even
    // since each encoded symbol occupies 2 bits
    mask[0] = 0x0000000000000000;
    // mask[0] = 0xF000000000000000 & MASK_0x1;;
    for (uint32_t i = 1; i < 64; i++)
    {
        mask[i] = mask[i-1] | (0x1UL << (64 - i));
    }
}

/* convert from encoded symbols to array of chars */
/* For example:
 * 10 00 11 11 00 01 00
 *  G  A  T  T  A  C  A */
void
decode_symbols(uint idx, char* seq_out, char* alphabet, uint bits_symbol, uint len)
{
  uint mask = (0x1LU << bits_symbol) - 1;

  for(uint i=0; i < len; i++)
      seq_out[i] = alphabet[(idx >> ((len-1-i)*bits_symbol)) & mask];

   seq_out[len] = 0;
}

int
generate_SFM(SFM_t *fmi, uint8_t** bwt, uint64_t len,
             char * alphabet, size_t alignment, uint64_t* end_char_pos)
{
  // Prologue generation
  fmi->len = len;
  fmi->alphabet = alphabet;
  fmi->end_char_pos = end_char_pos;
  // len+1 because ep+1 is used
  fmi->n_entries = ceil_uint_div(len+1, D_VAL);
  fmi->C = calloc(K2_SYMBOLS + 1, sizeof(uint64_t));
  if (fmi->C == NULL)
  {
    fprintf(stderr, "Error when malloc fm-index memory\n");
    return -1;
  }

  // Rocc memory allocation
  fmi->entries = aligned_alloc(alignment, fmi->n_entries*sizeof(SFM_entry_t)*K2_SYMBOLS);
  if (fmi->entries == NULL)
  {
    fprintf(stderr, "Error when malloc fm-index memory\n");
    return -1;
  }

  // ROcc initialization
  for(uint64_t i=0; i < fmi->n_entries*K2_SYMBOLS; i++)
  {
    fmi->entries[i].counter = 0;
    fmi->entries[i].data = 0;
  }

  // Rocc generation
  uint32_t entry_id, word_offset, entry_offset;
  char c0, c1, c;

  for(uint64_t i=0; i<len; i++)
  {
    entry_id     = i / D_VAL;
    entry_offset = i % D_VAL;
    word_offset  = entry_offset % 64;

    if (entry_offset == 0)  // if (i % D_VAL == 0)
    {
      for(int j=0; j < K2_SYMBOLS; j++)
        fmi->entries[entry_id*K2_SYMBOLS + j].counter = fmi->C[j+1];
    }

    /* BWT[0]: least significant bits, BWT[1]: most significant bits */
    /* c: BWT[1] BWT[0] */
    c0 = read_char_from_buffer(bwt[0], BITS_PER_SYMBOL, i*BITS_PER_SYMBOL);
    c1 = read_char_from_buffer(bwt[1], BITS_PER_SYMBOL, i*BITS_PER_SYMBOL);
    c  = (c1 << BITS_PER_SYMBOL) | c0;

    // Write symbol in bitmap
    if ((fmi->end_char_pos[0] != i) && (fmi->end_char_pos[1] != i))
    {
      /* word_offset = 0 -> MSB data[] */
      uint64_t mask = 0x1LU << (63 - word_offset);
      fmi->entries[entry_id*K2_SYMBOLS + c].data |= mask;
      // Update C table
      fmi->C[c+1]++;
    }
  }

  // end+1 correction
  if (len % D_VAL == 0)
  {
    for(uint32_t i=0; i < K2_SYMBOLS; i++)
    {
      fmi->entries[(fmi->n_entries-1)*K2_SYMBOLS + i].counter = fmi->entries[(fmi->n_entries-2)*K2_SYMBOLS + i].counter;
      fmi->entries[(fmi->n_entries-1)*K2_SYMBOLS + i].data    = fmi->entries[(fmi->n_entries-2)*K2_SYMBOLS + i].data;
    }
  }
  
  // $ Symbol adjustments at C table
  // $ c0
  fmi->C[0] = 1;

  // c1 $
  c1 = read_char_from_buffer(bwt[1], BITS_PER_SYMBOL, end_char_pos[0]*BITS_PER_SYMBOL);
  fmi->last_char = c1;
  fmi->C[c1*4]++;  // fmi->C[c1 << BITS_PER_SYMBOL]++; ;
  // printf("bwt1[end_char_pos]=%u\n", c1);

  // C Table accumulation
  for(uint32_t i=1; i <= K2_SYMBOLS; i++)
    fmi->C[i] += fmi->C[i-1];
    
  // Adding C to the Occ entries
  for(uint64_t i=0; i<fmi->n_entries; i++)
  {
    for(uint32_t j=0; j < K2_SYMBOLS;j++)
      fmi->entries[i*K2_SYMBOLS + j].counter += fmi->C[j];
  }
  
  // Encoding tables generation
  generate_SFM_encoding_table(fmi, &(fmi->encoding_table));
  generate_SFM_encoding_table2(fmi, &(fmi->encoding_table2));
  
  // LUT generation
  // first dimension of LUT[x][y] indexed with len % 3 (x = len % 3)
  generate_SFM_LUT(fmi, 6, &(fmi->LUT[0]));
  generate_SFM_LUT(fmi, 5, &(fmi->LUT[1]));
  
  return 0;
}

int dump_SFM(SFM_t *fmi)
{
  printf("** SFM **\n");

  // Prologue
  printf("- start: %s\n", fmi->start);
  printf("- len: %lu\n", fmi->len);
  printf("- alphabet: %s\n", fmi->alphabet);
  printf("- end_char_pos: %lu\n", fmi->end_char_pos[0]);
  
  // C array
  dump_C_SFM(fmi->C);

  // Occ entries

  printf("- ROcc entries: %lu\n", fmi->n_entries);
  printf("                "); 
  for(int i = 0; i < D_VAL; i++)
    printf("%2u ", i);
  printf("\n"); 

  for(uint32_t i = 0; i < fmi->n_entries; i++)
  {
    for(int j = 0; j < K2_SYMBOLS; j++)
    {
        printf("[%03u][%c%c] ", i*K2_SYMBOLS + j, fmi->alphabet[(j & 0xC) >> 2], fmi->alphabet[j & 0x3]);
        printf("%3lu ", fmi->entries[i*K2_SYMBOLS + j].counter);
        printf("| ");
        // printf("%2lx ", fmi->entries[i].data);
        for(int k = 0; k < D_VAL; k++)
        {
            printf("%2lu ", (fmi->entries[i*K2_SYMBOLS + j].data >> (63 - k)) & 0x1);
        }
        printf("\n"); 
    }
    printf("----------------------------------------------------------------------------------------------"); 
    printf("----------------------------------------------------------------------------------------------\n"); 
  }

  // 6-char LUT
  uint lut_entries = pow(SYMBOLS, 6);
  printf("- LUT entries: %u\n", lut_entries);
  for(uint i = 0; i < lut_entries; i++)
    printf("  [%03u] %u - %u\n", i, (fmi->LUT[0])[i].start, (fmi->LUT[0])[i].end);

  return 0;
}

int write_SFM(const char* file, SFM_t *fmi)
{
  FILE *f = fopen(file, "w");
  if (f == NULL)
  {
    fprintf(stderr, "Cannot open file %s \n", file);
    return -1;
  }

  // Write prologue
  fwrite(fmi->start, sizeof(char), 500, f);
  fwrite(&fmi->len, sizeof(fmi->len), 1, f);
  fwrite(fmi->alphabet, sizeof(char), SYMBOLS, f);
  fwrite(fmi->end_char_pos, sizeof(uint64_t), KSTEPS, f);

  // Write C array
  fwrite(fmi->C, sizeof(uint64_t), K2_SYMBOLS+1, f);

  // Write Occ entries
  fwrite(&fmi->n_entries, sizeof(fmi->n_entries), 1, f);
  fwrite(fmi->entries, sizeof(SFM_entry_t), K2_SYMBOLS*fmi->n_entries, f);

  // Write encoding table
  fwrite(fmi->encoding_table, sizeof(unsigned char), 256, f);
  fwrite(fmi->encoding_table2, sizeof(unsigned char), 256*256, f);

  // Write LUTs
  // LUT table [0], [1] indexed with len % 2
  uint lut_entries6 = pow(SYMBOLS, 6);
  fwrite(fmi->LUT[0], sizeof(LUT_entry_t), lut_entries6, f);

  uint lut_entries5 = pow(SYMBOLS, 5);
  fwrite(fmi->LUT[1], sizeof(LUT_entry_t), lut_entries5, f);

  fclose(f);
  return 0;
}


int load_SFM(const char* file, SFM_t * fmi)
{
  FILE *f = fopen(file, "rb");
  if (f == NULL)
  {
    fprintf(stderr, "Cannot open file %s\n", file);
    return -1;
  }
  long fr;

  // Read 500 chars (TTGGATCTATGCTTCTGGT...)
  fmi->start = malloc(sizeof(char)*501);
  fr = fread(fmi->start, sizeof(char), 500, f);
  if (fr != 500)
  {
    fprintf(stderr, "Error at fread(start)\n");
    return -1;
  }
  fmi->start[500]=0;

  // Read prologue
  fr = fread(&fmi->len, sizeof(fmi->len), 1, f);
  if (fr != 1)
  {
    fprintf(stderr, "Error at fread(len)\n");
  }

  fmi->C = calloc(K2_SYMBOLS+1, sizeof(uint64_t));
  fmi->alphabet = malloc(SYMBOLS*sizeof(char));
  if ((fmi->C == NULL) || (fmi->alphabet == NULL))
  {
    fprintf(stderr, "Error when malloc fm-index memory\n");
    exit(1);
  }
  fr = fread(fmi->alphabet, sizeof(char), SYMBOLS, f);
  if(fr != SYMBOLS)
  {
    fprintf(stderr, "Error at fread(alphabet)\n");
    exit(1);
  }
  fmi->end_char_pos = malloc(sizeof(uint64_t)*KSTEPS);
  fr = fread(fmi->end_char_pos, sizeof(uint64_t), KSTEPS, f);
  if (fr != KSTEPS)
  {
    fprintf(stderr, "Error at fread(end_char_pos)\n");
    exit(1);
  }

  //read C array
  fr = fread(fmi->C, sizeof(uint64_t), K2_SYMBOLS+1, f);
  if (fr != K2_SYMBOLS+1)
  {
    fprintf(stderr, "Error at fread(C)\n");
    exit(1);
  }

  // read Occ entries
  fr = fread(&fmi->n_entries, sizeof(fmi->n_entries), 1, f);
  if (fr != 1)
  {
    fprintf(stderr, "Error at fread(n_entries)\n");
    exit(1);
  }

  uint64_t bytes_SFM = K2_SYMBOLS*fmi->n_entries*sizeof(SFM_entry_t);

#ifdef KNL
  fr = hbw_posix_memalign_psize((void**)&fmi->entries, 64, bytes_SFM, HBW_PAGESIZE_1GB);
  if(fr != 0)
  {
    fprintf(stderr, "%li Error at hbwmalloc (Entries)\n", fr);
    exit(1);
  }
#elif defined (HUGEPAGES)
  fmi->entries = mmap(ADDR, bytes_SFM, PROTECTION, HUGE_PAGE_FLAGS, 0, 0);
  if (fmi->entries == MAP_FAILED)
  {
    fflush(stdout);   
    perror("  mmap");
    printf("  Failed to allocate SFM %.1f%ciB in huge pages.\n  Trying conventional memory ...",
           bytes_SFM/V_1MB > 1024.0? bytes_SFM/V_1GB : bytes_SFM/V_1MB,
           bytes_SFM/V_1MB > 1024.0? 'G' : 'M');
    fr = posix_memalign((void**)&fmi->entries, ALIGN_2MB, bytes_SFM);
    if (fr != 0)
    {
      fprintf(stderr, "%li Error at malloc (Entries)\n", fr);
      exit(1);
    }
    printf("OK\n");
  }
  else
  {
    printf("  %.1f%ciB SFM allocated in huge pages\n",
           bytes_SFM/V_1MB > 1024.0? bytes_SFM/V_1GB : bytes_SFM/V_1MB,
           bytes_SFM/V_1MB > 1024.0? 'G' : 'M');
  }
#else
  fr = posix_memalign((void**)&fmi->entries, ALIGN_2MB, bytes_SFM);
  if (fr != 0)
  {
    fprintf(stderr, "%li Error at malloc (Entries)\n", fr);
    exit(1);
  }
  // fprintf(stderr, "SFM NO mapeada en huge page\n"); 
#endif

  fr = fread(fmi->entries, sizeof(SFM_entry_t), fmi->n_entries*K2_SYMBOLS, f);
  if (fr != (long) fmi->n_entries*K2_SYMBOLS)
  {
    fprintf(stderr, "Error at fread (Entries)\n");
    exit(1);
  }

  // read encoding table  
  fmi->encoding_table = (unsigned char*)malloc(sizeof(unsigned char)*256);
  fr =fread(fmi->encoding_table, sizeof(unsigned char), 256, f);
  if (fr != 256)
  {
    fprintf(stderr, "Error at fread (Encoding table)\n");
    exit(1);
  }
  
  fmi->encoding_table2 = (unsigned char*)malloc(sizeof(unsigned char)*256*256);
  fr = fread(fmi->encoding_table2, sizeof(unsigned char), 256*256, f);
  if (fr != 256*256)
  {
    fprintf(stderr, "Error at fread (Encoding table)\n");
  }  

  // read LUTs
  uint lut_entries6 = pow(SYMBOLS, 6);
  fmi->LUT[0] = (LUT_entry_t *) malloc(lut_entries6*sizeof(LUT_entry_t));
  if (fmi->LUT[0] == NULL)
  {
    fprintf(stderr, "Error at LUT6 malloc. \n");
    exit(1);
  }
  fr = fread(fmi->LUT[0], sizeof(LUT_entry_t), lut_entries6, f);
  if (fr != lut_entries6)
  {
    fprintf(stderr, "Error at fread (LUT6)\n");
      exit(1);
  }
  uint lut_entries5 = pow(SYMBOLS, 5);
  fmi->LUT[1] = (LUT_entry_t *) malloc(lut_entries5*sizeof(LUT_entry_t));
  if (fmi->LUT[1] == NULL)
  {
    fprintf(stderr, "Error at LUT5 malloc. \n");
    exit(1);
  }
  fr = fread(fmi->LUT[1], sizeof(LUT_entry_t), lut_entries5, f);
  if (fr != lut_entries5)
  {
    fprintf(stderr, "Error at fread (LUT5)\n");
    exit(1);
  }
  fclose(f);
  return 0;
}

int
generate_SFM_encoding_table(SFM_t *fmi, uint8_t** out)
{
    *out = (uint8_t*) malloc(256*sizeof(uint8_t));

    for(uint i=0; i < 256; i++)
      (*out)[i] = -1;

    for(uint i=0; i < SYMBOLS; i++)
      (*out)[(uint8_t) fmi->alphabet[i]] = i;

    return 0;
}

int
count_SFM(SFM_t *fmi, const char* seq, uint len, uint64_t * start, uint64_t * end)
{
    uint step_len = ((len-1)/KSTEPS);
    uint8_t ch1, ch2, ch;

    if (len % 2 == 0) // If the length of the sequence is even
    { // Search the last two characters
      ch1 = fmi->encoding_table[(uint)seq[len-1]];
      ch2 = fmi->encoding_table[(uint)seq[len-2]];
      ch = ch1 + (ch2 << BITS_PER_SYMBOL);
      *start = fmi->C[(uint)ch];
      *end   = fmi->C[ch+1];
    }
    else
    { // Search the last character individually
      ch = fmi->encoding_table[(uint)seq[len-1]];
      *start = fmi->C[ch*4];
      *end   = fmi->C[(ch+1)*4];
      if (ch == fmi->last_char)
        (*start)--;
    }

    for(int i=step_len-1; i>=0; i--)
    {
      // Encode the char to search
      ch1 = fmi->encoding_table[(uint)seq[2*i + 1]];
      ch2 = fmi->encoding_table[(uint)seq[2*i + 0]];
      ch = ch1 + (ch2 << BITS_PER_SYMBOL);

      // Start and end interval update
      *start = k2_LF(ch, *start);
      *end   = k2_LF(ch, *end);
    }
    (*end)--;
    return 0;
}

int
generate_SFM_encoding_table2(SFM_t *fmi, uint8_t** out)
{
    uint symbol;
    char encoded_symbol, ch1, ch2;

    *out = (uint8_t*) malloc(256*256*sizeof(uint8_t));
    for(uint i = 0; i < 256*256; i++)
      (*out)[i] = -1;

    for(uint i = 0; i < SYMBOLS; i++)
    {
      for(uint j = 0; j < SYMBOLS; j++)
      {
        ch1 = fmi->alphabet[j];
        ch2 = fmi->alphabet[i];
        symbol = ch1 + (ch2 << 8);
        encoded_symbol = i + (j << BITS_PER_SYMBOL);
        (*out)[symbol] = encoded_symbol;
      }
    }
    return 0;
}

int
generate_SFM_LUT(SFM_t * fmi, uint n_chars, LUT_entry_t ** LUT)
{
  uint lut_entries = pow(SYMBOLS, n_chars);
  uint64_t start, end;
  char * seq = malloc((n_chars+1)*sizeof(char));
  *LUT = (LUT_entry_t *) malloc(lut_entries*sizeof(LUT_entry_t));

  // Data initializing
  pfmi = fmi;
  ROcc = fmi->entries;
  mask_init(mask_64b);

  if ((seq == NULL) || (*LUT == NULL))
  {
    fprintf(stderr, "Error at LUT malloc.\n");
    return -1;
  }

  for(uint i=0; i < lut_entries; i++)
  {
    decode_symbols(i, seq, fmi->alphabet, BITS_PER_SYMBOL, n_chars);
    count_SFM(fmi, seq, n_chars, &start, &end);
    (*LUT)[i].start = start;
    (*LUT)[i].end = end;
  }
  free(seq);
  return 0;
}

void
free_SFM(SFM_t * fmi)
{
  free(fmi->encoding_table);
#ifdef KNL
  hbw_free(fmi->entries);
#elif defined (HUGEPAGES)
  munmap(fmi->entries, (long)fmi->n_entries*K2_SYMBOLS*sizeof(SFM_entry_t));
#else
  free(fmi->entries);
#endif
  free(fmi->alphabet);
  free(fmi->start);
  free(fmi->end_char_pos);
  free(fmi->C);
  free(fmi->LUT[0]);
  free(fmi->LUT[1]);
}
