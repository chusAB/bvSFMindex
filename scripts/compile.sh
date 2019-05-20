#!/bin/bash

# default values
version="k2d64bv"
arch=knl
comp=gcc
nseqs=4
pfetch=sp

perf=0
build=0
iaca=0

arch_flag=0
comp_flag=0

while getopts "a:v:c:s:p:w:bih" opt; do
  case $opt in
    a) 
      # echo "architecture -> $OPTARG"
      arch=$OPTARG
      ;;
    v) 
      # echo "version -> $OPTARG"
      version=$OPTARG
      ;;
    c) 
      # echo "compiler -> $OPTARG"
      comp=$OPTARG
      ;;
    s) 
      # echo "number of overlapping sequences -> $OPTARG"
      nseqs=$OPTARG
      ;;
    p)
      # echo "prefetch type (sp|dp) -> $OPTARG"
      pfetch=$OPTARG
      ;;
    w)
      # echo "hardware counters with perf -> $OPTARG"
      perf=$OPTARG
      ;;
    b) 
      # echo "build compilation -> $OPTARG"
      build=1
      ;;
    i) 
      # echo "compilation with IACA support -> $OPTARG"
      iaca=1
      ;;
    h)
      echo "use:"
      echo "$0 -v version  -c compiler (icc,gcc-6,gcc-7) -a architecture (knl,bdw,skx,ivb,native)"
      echo "ejemplo:"
      echo "$0 -v k2d64bv -c gcc-7  -a knl"
      exit
      ;;
    \?)
      echo "invalid option: -$OPTARG"
      ;;
    :)
      echo "-$OPTARG option requires a parameter"
      exit 1
      ;;
  esac
done

printf "compiling ${version}.${arch}.${comp}.${nseqs}seq.${pfetch} ... "
#echo "compiling ${version} with ${comp} for ${arch} architecture..."

case ${arch} in
    gen) arch_flag=0
          ;;
    nat) arch_flag=1
          ;;
    knl) arch_flag=2
         ;;
    ivb) arch_flag=3
          ;;
    hsw) arch_flag=4
          ;;
    bdw) arch_flag=5
          ;;
    skl) arch_flag=6
          ;;
    skx) arch_flag=7
          ;;
    *)   arch_flag=0
esac

case ${comp} in
    icc)   comp_flag=0
           ;;
    gcc)   comp_flag=1
           ;;
    gcc-6) comp_flag=2
           ;;
    gcc-7) comp_flag=3
           ;;
    *)     comp_flag=1
esac

cd ../${version}
make  a=${arch_flag} c=${comp_flag} s=${nseqs} p=${pfetch} clean
make  a=${arch_flag} c=${comp_flag} s=${nseqs} p=${pfetch} w=${perf} i=${iaca} ${version}_fcount.${arch}.${comp}.${nseqs}seq.${pfetch} 

if [ ${build} -eq 1 ]; then
    make  a=${arch_flag} c=${comp_flag} ${version}_build.${arch}.${comp}
fi

printf "OK\n"

