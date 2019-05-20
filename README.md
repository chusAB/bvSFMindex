# bvSFM

[bvSFM] (bit-vector sampled FM-index) is a tool for sequence alignment. Specifically, it implements an exact search
algorithm that counts the number of matches of arbitrary length reads on a reference genome.
[bvSFM] indexes a genome with an [FM Index][FM Index Wiki] (based on the [Burrows-Wheeler Transform] or [BWT]).
[FM Index] is a compact data structure suitable for fast matches of short reads to large reference genomes.
For the human genome, its memory footprint is typically around 3.2 gigabytes of RAM.

[bvSFM] uses an optimized FM-index data structure layout and codification that
packs all relevant data needed in a query step within a single cache block,
minimizing the memory bandwidth demand.
[bvSFM] achieves best results when executed on multicore systems integrating high bandwidth memory,
for instance an Intel Xeon Phi processor [KNLa][KNLb] (codenamed Knights Landing, or KNL).

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


## Getting Started

[bvSFM] is available from this [GitHub repository].

Build the bvSFM tools by running GNU `make` with no arguments:

    $ cd k2d64bv
    $ make


## Prerequisites

### Required libraries

1) The [memkind] library is used to allocate the bvSFM-index in the KNL High Bandwidth Memory.
   Namely, the hbwmalloc API is used (hbw_posix_memalign_psize() function).

        $ mkdir memkindtmp
        $ cd memkindtmp
        $ git clone https://github.com/memkind/memkind.git
        $ cd memkind
        $ ./configure
        $ make
        $ sudo make install


2) libdivsufsort is a software library that implements a lightweight suffix array construction algorithm.
   We use this library to compute the [Burrows-Wheeler Transform] or [BWT] needed to build the bvSFM index.
   Hence, it is required to compile the `bvSFM` indexer.
   It can be installed in several ways:

    a) Directly from the [divsufsort library] source code:

        $ git clone https://github.com/y-256/libdivsufsort.git
        $ cd libdivsufsort
        $ mkdir build
        $ cd build
        $ cmake -DCMAKE_BUILD_TYPE="Release" -DBUILD_DIVSUFSORT64:BOOL=ON -DCMAKE_INSTALL_PREFIX="/usr/local" ..
        $ make
        $ sudo make install

    b) With the [Succinct Data Structure Library (SDSL)]. To build and install SDSL, follow the instructions at its `README` file:

        $ mkdir sdsltmp
        $ cd sdsltmp
        $ git clone https://github.com/simongog/sdsl-lite
        $ cd sdsl-lite/
        $ sudo ./install.sh /usr/local/
        [...]
        -- Installing: /usr/local/lib/libsdsl.a
        SUCCESS: sdsl was installed successfully!
        The sdsl include files are located in '/usr/local/include'.
        The library files are located in '/usr/local/lib'.
         
        Sample programs can be found in the examples-directory.
        A program 'example.cpp' can be compiled with the command: 
        g++ -std=c++11 -DNDEBUG -O3 [-msse4.2] \
           -I/usr/local/include -L/usr/local/lib \
           example.cpp -lsdsl -ldivsufsort -ldivsufsort64
         
        Tests in the test-directory
        A cheat sheet in the extras/cheatsheet-directory.
        Have fun!


### Optional libraries


1) If you want the application to print system information (number of nodes, number of cores)
   process information (node and core where the application is running),
   and huge page status, libnuma is required.

   To install libnuma in an CentOS/Red Hat system:

        $ sudo yum install numactl numactl-devel

   In an Debian/Ubuntu system:
   
        $ sudo apt-get install libnuma-dev 

   In case libnuma is available, compile with the following Makefile argument:

        $ make n=1


2) If you want to collect hardware counters when executing the search algorithm,
   the perf monitoring tools are required. To install perf in a Debian/Ubuntu system: 

        $ apt-get install linux-tools-common linux-tools-generic linux-tools-`uname -r`
   
   In case perf is available, compile with the following Makefile argument:

        $ make p=1


3) If you want to analyze with IACA the main loop of the search algorithm that
   counts the number of matches, compile with the following Makefile argument:

        $ make i=1
        
   IACA is supposed to be installed in the /opt/intel/iaca-lin64/ directory.
   It is hardcoded in the k2d64bv/k2d64bv.c file.


### Installing

Clone the repository:

    $ git clone https://github.com/memkind/memkind.git


To compile programs, just run `make` from the version directory.

    $ cd k2d64bv
    $ make

See the `MANUAL.md` file for more details of the compilation options.


### Compiling and running scripts

Several scripts are provided to automate the process of executing several configurations.
Edit the script `script/run_all.sh` to specify the configurations you want to execute.


## Developing new bvSFM versions

Create a new directory, for instance, `k2d96bv`, and create the corresponding
 `Makefile`, `k2d96bv.c`, `k2d96bv.h`, `k2d96bv_build.c`, and `k2d96bv_fcount.c` files.


## Authors

* **José-Manuel Herruzo-Ruiz** - jmherruzo@uma.es - *Initial work* - [jmherruzo web page]
* **Jesús Alastruey-Benedé** - jalastru@unizar.es - [chus web page]
* **Pablo Ibáñez-Marín** - imarin@unizar.es - [imarin web page]


## License

This project is licensed under the GNU General Public License v3.0 - see the `COPYING` file for details.



[BWT]:                                                http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[bvSFM]:                                              https://ieeexplore.ieee.org/document/8566001
[Burrows-Wheeler Transform]:                          http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[chus web page]:                                      http://webdiis.unizar.es/u/chus/
[FM Index]:                                           http://portal.acm.org/citation.cfm?id=796543
[FM Index Wiki]:                                      http://en.wikipedia.org/wiki/FM-index
[GitHub repository]:                                  https://github.com/chusAB/bvSFMindex
[GPLv3 license]:                                      http://www.gnu.org/licenses/gpl-3.0.html
[imarin web page]:                                    http://webdiis.unizar.es/u/imarin/
[jmherruzo web page]:                                 https://github.com/jmherruzo
[KNLa]:                                               https://ieeexplore.ieee.org/document/7453080
[KNLb]:                                               https://www.elsevier.com/books/intel-xeon-phi-processor-high-performance-programming/jeffers/978-0-12-809194-4
[k-step FM-Index]:                                    https://github.com/achacond/k-step_FM-index
[divsufsort library]:                                 https://github.com/y-256/libdivsufsort/
[memkind]:                                            https://github.com/memkind/memkind
[Succinct Data Structure Library (SDSL)]:             https://github.com/simongog/sdsl-lite
