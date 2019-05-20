#!/bin/bash

# uso:
#    ./run_all.sh

version_list=("k2d64bv")
# version_list=("k2d96bv")

arch_list=("gen")
# arch_list=("gen" "nat" "knl" "ivb" "hsw" "bdw" "skl" "skx")

comp_list=("gcc")
# comp_list=("icc" "gcc-6" "gcc-7" )

# number of hardware threads
nthreads_list=("1" "2" "3" "4")
# nthreads_list=("1" "2" "4" "8" "16" "28" "56")

# number of overlapped searches
# skl
nseq_list=("16")
# nseq_list=("16" "24" "32" "40")
# knl
# nseq_list=("2" "4" "6" "8")

# types of software prefetch
pfetch_list=("sp")
# pfetch_list=("sp" "dp")

# genome reference
gref=lambda_virus

# file containing the sequences to search
file_seq=reads_1.fasta

# collect hw counters with perf
perf=0

# loop over the versions
for version in "${version_list[@]}"
do
    # loop over the architectures
    for arch in "${arch_list[@]}"
    do
        # loop over the compilers
        for comp in "${comp_list[@]}"
        do
            # loop over the overlapped searches
            for nseq in "${nseq_list[@]}"
            do
                # loop over the prefetches
                for pfetch in "${pfetch_list[@]}"
                do
                    # -b: ejecucion de build, -d: ejecucion de count
                    ./compile.sh -v $version -a $arch -c $comp -s $nseq -p $pfetch -w $perf
                    # loop over the threads
                    for nthreads in "${nthreads_list[@]}"
                    do
                        ./run.sh  -v $version -a $arch -c $comp -t $nthreads -s $nseq -p $pfetch -w $perf -g $gref -f $file_seq
                    done
                done
            done
        done
    done
done
