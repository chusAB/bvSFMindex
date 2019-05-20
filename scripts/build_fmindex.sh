#!/bin/bash

# valores por defecto
version="k2d64bv"
arch=nat
comp=icc
gref=GRCh38_1MB

while getopts "a:v:c:g:x:h" opt; do
  case $opt in
    a) 
      # echo "especificada architectura -> $OPTARG"
      arch=$OPTARG
      ;;
    v) 
      # echo "especificada version -> $OPTARG"
      version=$OPTARG
      ;;
    c) 
      # echo "especificado compilador -> $OPTARG"
      comp=$OPTARG
      ;;
    g) 
      # echo "especificada referencia -> $OPTARG"
      gref=$OPTARG
      ;;
    x)
      wt=$OPTARG
      ;;
    h)
      echo "use:"
      echo "$0 -v version  -c compiler -a architecture -g genome_reference"
      echo "example:"
      echo "$0 -v k2d64bv -c gcc-7 -a knl -g GRCh38"
      exit
      ;;
    \?)
      echo "invalid optiom: -$OPTARG"
      ;;
    :)
      echo "-$OPTARG option requires a parameter"
      exit 1
      ;;
  esac
done

echo "generating index ${version} for reference ${gref} with version compiled with ${comp} for the ${arch} architecture ..."
PREFIX=../bin
echo "executing ${PREFIX}/${version}_build.${arch}.${comp}  ../references/${gref} ... "
${PREFIX}/${version}_build.${arch}.${comp}  ../references/${gref}
