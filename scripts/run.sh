#!/bin/bash

# default values
version="k2d64bv"
arch=knl
comp=gcc
nthreads=256
nseqs=4
pfetch=sp

perf=0
outdir=runs

gref=GRCh38
file_seq=seq_20M.fa

while getopts "a:v:c:t:s:p:w:g:f:h" opt; do
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
    t) 
      # echo "number of threads -> $OPTARG"
      nthreads=$OPTARG
      ;;
    s) 
      # echo "number of overlapped sequences -> $OPTARG"
      nseqs=$OPTARG
      ;;
    p)
      # echo "prefetch type (sp|dp) -> $OPTARG"
      pfetch=$OPTARG
      ;;
    g) 
      # echo "reference genome -> $OPTARG"
      gref=$OPTARG
      ;;
    f) 
      # echo "file with sequences to search -> $OPTARG"
      file_seq=$OPTARG
      ;;
    w)
      # echo "collect hardware counters with perf -> $OPTARG"
      perf=$OPTARG
      ;;
    h)
      echo "uso:"
      echo "$0 -v version  -c compiler -a architecture -g reference_genome -f sequence"
      echo "example:"
      echo "$0 -v k2d64bv -c gcc-7 -a knl -g GRCh38_1GB -f seq_10M.fa"
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

if [ ${perf} -eq 1 ]; then
    outdir=perf
fi

mkdir -p ../${outdir}
fseq=${file_seq%.*}
nth=`printf "%02d\n" $nthreads`
nsq=`printf "%02d\n" $nseqs`

outfile=../${outdir}/${version}.${arch}.${comp}.${nth}th.${nsq}seq.${pfetch}.${gref}.${fseq}.`date +%Y%m%d.%02H%M%S`.`hostname | cut -f1 -d.`.txt
touch ${outfile}
echo -n "Test machine: "  >> ${outfile}
hostname >> ${outfile}
echo -n "Model name:"  >> ${outfile}
cat /proc/cpuinfo | grep 'model name' | uniq | cut -f2 --delimiter=: >> ${outfile}
#cat /proc/cpuinfo | grep 'model name' | uniq | cut -c 13-40 >> ${outfile} 
#TODO: dump number of cores
echo -n "Date: " >> ${outfile} 
LANG=en_EN date >> ${outfile}

PREFIX=../bin
bin=${PREFIX}/${version}_fcount.${arch}.${comp}.${nseqs}seq.${pfetch}
reffile=../references/${gref}.${version}.fmi
seqfile=../sequences/${file_seq}
echo -n "executing ${bin} -f ${reffile}  -s ${seqfile} -t ${nthreads} ... "
${bin} -f ${reffile} -s ${seqfile} -t ${nthreads} >> ${outfile} 2>&1
# If the aligner is executed in a multiprocessor system,
# best results are obtained if all the threads are executed in the same processor
# For instance, for a 2xIntel Xeon Gold 5120 system:
# taskset -c  0-13,28-41 ${bin} -f ${reffile} -s ${seqfile} -t ${nthreads} >> ${outfile} 2>&1
# or
# taskset -c 14-27,42-55 ${bin} -f ${reffile} -s ${seqfile} -t ${nthreads} >> ${outfile} 2>&1

# OMP_PROC_BIND=close - bind threads close to the master thread while still distributing threads for load balancing.
# OMP_PLACES='sockets(1)' - only allow the application to run on the cores provided by a single CPU socket
# OMP_PROC_BIND=close OMP_PLACES='sockets(1)'  ${bin} ../References/${gref}.${version}.fmi ../Sequences/${file_seq} >> ${outfile} 2>&1
# KMP_AFFINITY=granularity=thread,compact ${bin} ${ref} ../Sequences/${file_seq} >> ${outfile} 2>&1
# KMP_AFFINITY=verbose,granularity=thread,compact ${bin} ../References/${gref}.${version}.fmi ../Sequences/${file_seq} >> ${outfile} 2>&1
# KMP_AFFINITY=verbose,granularity=core,compact ${bin} ../References/${gref}.${version}.fmi ../Sequences/${file_seq} >> ${outfile} 2>&1
# taskset -c 14-27,42-55 ${bin} ../References/${gref}.${version}.fmi ../Sequences/${file_seq} >> ${outfile} 2>&1
echo >> ${outfile}
printf "OK\n\n"
