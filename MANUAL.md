

# Introduction

[bvSFM] (bit-vector sampled FM-index) is a tool for sequence alignment. Specifically, it implements an exact search
algorithm that counts the number of matches of arbitrary length reads on a reference genome.
[bvSFM] indexes a genome with an [FM Index][FM Index Wiki] (based on the [Burrows-Wheeler Transform] or [BWT]).
[FM Index] is a compact data structure suitable for fast matches of short reads to large reference genomes.
For the human genome, its memory footprint is typically around 3.2 gigabytes of RAM.
The sequences to align have to be stored in a [FASTA] file.

[bvSFM] uses an optimized FM-index data structure layout and codification that
packs all relevant data needed in a query step within a single cache block,
minimizing the memory bandwidth demand.
[bvSFM] achieves best results when executed on multicore systems integrating high bandwidth memory,
for instance an Intel Xeon Phi processor [KNLa][KNLb] (codenamed KnightsLanding, or KNL).

[bvSFM] is distributed under the [GPLv3 license], and it runs on the command line under Linux.

If you use [bvSFM] for your published research, please cite the following paper:

J.M. Herruzo, S. González, P. Ibáñez, V. Viñals, J. Alastruey-Benedé, and Óscar Plata.
Accelerating Sequence Alignments Based on FM-Index Using the Intel KNL Processor.
IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB 2019).
DOI: 10.1109/TCBB.2018.2884701 

[bvSFM] can be used as a CPU benchmark. It performs intensive integer computation.
It can also be used to measure random memory access bandwidth.


## Exact Matching Using FM-Index

Given a pattern Q[1..p], the FM-index allows to find all occurrences of Q in a reference text T.
The search process takes two steps: Count and Locate. The first step is a rank query
process that calculates the number of occurrences of Q in T.
The second step uses the indexes of these rows to access the suffix array, where it
finds the position of every occurrence of Q in the text T.
The suffix array is usually a very large compressed data structure.
However, its size for the human genome is about 12 GB (3 gigabases x 4 bytes),
so it can be stored without compression in modern systems.
This way, the Locate step is very simple, it only requires an access to the suffix array.
[bvSFM] only implements the Count step.


## k2d64bv version

The Sampled FM-index (SFM) stores one set of counters out of every d sets in the original Occ structure.
This variant introduces a trade-off between memory footprint and computing cost.

The [k-step FM-Index] searches k symbols in a query in a single step.
This version reduces computing cost and improves data locality of the sampled version
but increases slightly memory footprint.

This repository releases a k2d64bv version, that is, a sampling factor (d) of 64
and a k-step value of 2.


# Obtaining bvSFM

[bvSFM] is available from this [GitHub repository].


## Building from source

Building [bvSFM] from source requires a GNU-like environment with GCC, GNU Make
and other basics. It should be possible to build [bvSFM] on most vanilla Linux
installations or on a Mac installation with [Xcode] installed.

First, clone the repository from the [GitHub repository].

    $ git clone git@github.com:chusAB/bvSFMindex.git

Then build the bvSFM tools by running GNU `make` with no arguments:

    $ cd k2d64bv
    $ make


The `Makefile` accepts several parameters. By default, it uses the following options:

- Compiler (c): `gcc`
- Number of overlapped sequence searches (s): 4
- Target architecture (a): generic
- Prefetch (p): dual (L1+L2)
- Huge pages support : yes
- Perf analysis (w): no
- IACA analysis (i): no

The default values can be overriden from the command line, for instance:

    $ make c=0 a=1

will the compile with `icc` for the native architecture.


[bvSFM] can be run on many threads by using OpenMP.


## Adding to PATH

By adding your new [bvSFM] directory to your [PATH environment variable], you
ensure that whenever you run a [bvSFM] tool from the command line, you will
get the version you just installed without having to specify the entire path.
This is recommended for most users. To do this, follow your operating system's
instructions for adding the directory to your [PATH].

If you would like to install [bvSFM] by copying the [bvSFM] executable files
to an existing directory in your [PATH], make sure that you copy all the
executables.


# The `bvSFM` aligner

`bvSFM` takes a bvSFM FM-index and a sequencing read file and outputs
the number of exact matches.

"Alignment" is the process by which we discover where the read sequences
can be found in the reference sequence. An "alignment" is a result from this
process, specifically: an exact alignment is a way of "lining up" all of the
characters in the read with some characters from the reference. For example:

      Read:             CTGCGAT
                        |||||||
      Reference: ACGACACCTGCGATCTCGACTCG

Where vertical bars show where aligned characters match.


## Alignment summary

When [bvSFM] finishes running, it prints messages summarizing what happened.
These messages are printed to the "standard output" ("stdout") file handle.
The summary might look like this:

    Occurrences found: 359
    Total LFOP: 2.00G
    Total time: 15.032853
    Raw throughput:  0.133 GLFOPS


## Command Line

### Usage

    ./k2d64bv_fcount  -f fmindex  -s sequences -t nthreads [-r runs]

         -f, --fmindex
             file storing the fm-index
         -s, --sequences
             file storing the sequences to search
         -t, --nthreads
             number of threads
         -r, --runs
             number of runs
         -h, --help
             show program usage


# The `bvSFM` indexer
===========================

`bvSFM_build` builds a bvSFM index from a genome reference. `bvSFM_build` outputs
a file with suffix `.fmi`. These file is needed to align reads to that reference.
The original reference file is no longer used by bvSFM once the index is built.

The bvSFM index is based on the [FM Index][FM Index Paper] of Ferragina and Manzini,
which in turn is based on the [Burrows-Wheeler] transform (BWT).
The algorithm that computes the BWT uses the [Succinct Data Structure Library (SDSL)]
to perform the suffix array calculation.


## Command Line

Usage:

    ./k2d64bv_build reference_file

The reference input file is a text file containing the genome reference.


# Getting started with bvSFM: Lambda phage example

bvSFM comes with some example files to get you started. The example files
are not scientifically significant; we use the [Lambda phage] reference genome
simply because it's short, and the reads were generated by a computer program,
not a sequencer. However, these files will let you start running bvSFM and
downstream tools right away.


## Indexing a reference genome

To create an index for the [Lambda phage] reference genome included with bvSFM,
go to the bvSFM directory and run:

    $ bin/bvSFM_build.version_suffix references/lambda_virus

The command should print many lines of output then quit.
When the command completes,  the current directory will contain a new file
containing the fm-index: `lambda_virus.k2d64bv.fmi`.

For instance:

    $ bin/k2d64bv_build.nat.gcc-7 references/lambda_virus 
    Reading FM-index file references/lambda_virus... OK
    Total time: 0.000s
    Getting BWT of 48502 characters... OK
    BWT Generated. Length: 48503
    Total time: 0.004s
    Encoding BWT... OK
    Total time: 0.001s
    ---Encoded C Array---(48503 elements)
    C[0]:  12334  11362  12820  11986 
    C[1]:  12334  11362  12820  11986 
    Compressing BWT... OK
     -> Compression ratio: 2.00 = 48503 / 24252 (original/compressed bytes)
    Total time: 0.000s
    Generating FM-index... OK -> FM-index written to file references/lambda_virus.k2d64bv.fmi
    FM-index time: 0.001s


## Aligning example reads

    $ bin/bvSFM_fcount.version_suffix -f references/lambda_virus.k2d64bv.fmi -s sequences/reads

This aligns a set of reads to the [Lambda phage] reference genome using
the index generated in the previous step. A summary is written to the console.
Actually, the summary is written to the "standard output" or "stdout" filehandle,
which is typically printed to the console.

For instance:

    $ bin/k2d64bv_fcount.nat.gcc-7.4seq.dp -f references/lambda_virus.k2d64bv.fmi -s sequences/reads_1.fasta -t 1
    Program version: 20190206
    -------------------------------------------------------------
    System info: 1 node, 4 cores (4 cores/node)
    Process info: running in node 0, core 2.
    -------------------------------------------------------------
    Loading FM-index...
      mmap: Cannot allocate memory
      Failed to allocate SFM 0.2MiB in huge pages.
      Trying conventional memory ...OK
    OK. Index loaded in 0.000485s
    Free hugepages   2MB   1GB
    node 0:            0    NA
    -------------------------------------------------------------
    Loading sequences...
    OK. 0.01 Msequences loaded in 0.004547s (2.199 Mseq/s)
    -------------------------------------------------------------
    Parameters
    - FM-index file: references/lambda_virus.k2d64bv.fmi
    - Index size: 0.0GiB (48503 characters)
    - Number of threads: 1
    - Overlapped sequences: 4
    - Sequence file: sequences/reads_1.fasta
    - Number of bases: 0.00 Gbases (1062398)
    - Number of sequences: 0.01 Mseq (10000)
    - Average sequence length: 106.24 bases
    - Number of sequences processed per thread: 0.0M (10000)
    -------------------------------------------------------------
    The kernel will be executed 5 times.
    The *best* throughput will be reported (excluding the first iteration).
    -------------------------------------------------------------
    Starting search... 
    OK
    Occurrences found: 1192 (in each run)
    Total LFOP: 0.01G (expected 0.01G)
    Total time: 0.022235
    Raw throughput:  0.478 GLFOPS
    -------------------------------------------------------------
    Best throughput:  0.528 GLFOPS (0.528)
    Avg. throughput:  0.524 GLFOPS (0.524)
    Best bandwidth:   9.235 GB/s (9.244)
    Avg. bandwidth:   9.170 GB/s (9.244)
    Throughputs: [0.354] 0.517 0.528 0.525 0.526
    Min time:  0.004 s (run #2)
    Avg time:  0.004 s
    Max time:  0.004 s (run #1)
    Times: [0.006] 0.004 0.004 0.004 0.004
    Standard deviation = 0.772%
    -------------------------------------------------------------


## Performance tuning

1.  Huge pages

    When sequences are searched in an index of the order of gigabytes,
    huge page may improve performance by reducing the number of TLB misses.
    Check the `huge_pages.md` file to see how to enable huge page support.


2.  Overlapped queries

    The exact matching algorithm is typically query latency bound. 
    However, the memory latency can be hidden by issuing a given number of different independent queries,
    that is, by overlapping the memory accesses of several queries.
    
    Current architectures support simultaneous multithreading (SMT), a technique
    that allows a single core to execute several interleaved independent execution flows (hardware threads).
    In this situation, the queries can be distributed among the hardware threads of a core.
    
    The recommended overlapping factor is 4 for KNL, and 20 for Broadwell and Skylake processors.


3.  Multicore/Multiprocessor
 
    If your system has multiple processors/cores, use `-t`

    The [`-t`] option causes [bvSFM] to launch a specified number of parallel
    search threads. Each thread runs on a different processor/core and all
    threads find alignments in parallel, increasing alignment throughput.

    If the aligner is executed in a multiprocessor system,
    best results are obtained if all the threads are executed in the same processor
    For instance, for a 2xIntel Xeon Gold 5120 system:

        taskset -c  0-13,28-41 bin/k2d64bv_fcount.nat.gcc-7.4seq.dp -f references/lambda_virus.k2d64bv.fmi -s sequences/reads_1.fasta -t 1

    Or:
    
        taskset -c 14-27,42-55 bin/k2d64bv_fcount.nat.gcc-7.4seq.dp -f references/lambda_virus.k2d64bv.fmi -s sequences/reads_1.fasta -t 1



# Acknowledgements

This document is based on the [BOWTIE2] manual.



[BWT]:                                                http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[bvSFM]:                                              https://ieeexplore.ieee.org/document/8566001
[bvSFM paper]:                                        https://ieeexplore.ieee.org/document/8566001
[BOWTIE2]:                                            https://github.com/BenLangmead/bowtie2/
[Burrows-Wheeler Transform]:                          http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[Burrows-Wheeler]:                                    http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[Download]:                                           http://webdiis.unizar.es/~chus/
[FASTA]:                                              https://en.wikipedia.org/wiki/FASTA
[FM Index Paper]:                                     http://portal.acm.org/citation.cfm?id=796543
[FM Index Wiki]:                                      http://en.wikipedia.org/wiki/FM-index
[GitHub repository]:                                  https://github.com/chusAB/bvSFMindex
[GPLv3 license]:                                      http://www.gnu.org/licenses/gpl-3.0.html
[KNLa]:                                               https://ieeexplore.ieee.org/document/7453080
[KNLb]:                                               https://www.elsevier.com/books/intel-xeon-phi-processor-high-performance-programming/jeffers/978-0-12-809194-4
[k-step FM-Index]:                                    https://github.com/achacond/k-step_FM-index
[Lambda phage]:                                       http://en.wikipedia.org/wiki/Lambda_phage
[PATH environment variable]:                          http://en.wikipedia.org/wiki/PATH_(variable)
[PATH]:                                               http://en.wikipedia.org/wiki/PATH_(variable)
[Performance tuning]:                                 #performance-tuning
[Xcode]:                                              http://developer.apple.com/xcode/
[manual section on index building]:                   #the-bvSFM-build-indexer
[multiseed heuristic]:                                #multiseed-heuristic
[obtain bvSFM]:                                       #obtaining-bvSFM
[sourceforge site]:                                   https://sourceforge.net/projects/bowtie-bio/files/bowtie2/
[Succinct Data Structure Library (SDSL)]:             https://github.com/simongog/sdsl-lite
