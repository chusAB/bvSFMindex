#!/bin/bash

# uso:
#    ./compile_all.sh

version_list=("k2d64bv")
# version_list=("k2d96bv")

arch_list=("gen")
# arch_list=("knl" "bdw" "skx" "ivb" "native")

comp_list=("gcc")
# comp_list=("icc" "gcc-6" "gcc-7" )

# number of overlapped searches
# skl
nseq_list=("16")
# nseq_list=("16" "24" "32" "40")
# knl
# nseq_list=("2" "4" "6" "8")

# types of software prefetch
pfetch_list=("sp")
# pfetch_list=("sp" "dp")

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
                    # -b: compile FMindex build program
                    # -d: compile FMindex count program
                    ./compile.sh -v $version -a $arch -c $comp -s $nseq -p $pfetch -w $perf
                done
            done
        done
    done
done
