#!/bin/bash

# use:
#    ./build_all.sh

version_list=("k2d64bv")
arch_list=("gen")
# arch_list=("knl" "bdw" "skx" "ivb" "native")

comp_list=("gcc" )
# comp_list=("icc" "gcc-6" "gcc-7" )

# reference genome
gref=lambda_virus
# gref=GRCh38

# loop over the versions
for version in "${version_list[@]}"
do
    # loop over the architectures
    for arch in "${arch_list[@]}"
    do
        # loop over the compilers
        for comp in "${comp_list[@]}"
        do
            ./compile.sh -v $version -a $arch -c $comp -b
            ./build_fmindex.sh -v $version -a $arch -c $comp -g $gref
            # echo -v $version -a $arch -c $comp
        done
    done
done
