#!/bin/bash

if [ $# != 2 ]; then
 echo "Usage: ./make_2adic.sh <range> <src_dir>"
 exit 1
else
    r=$1
    d=$2
fi

echo 'running Magma 2-adic script on '${r} in ${d}

# Here we need the label in 3 parts then the a-invariant list, space-separated
# - in curvedata files this is fields 10,2,3,11
# - in allcurves files (as output by process_raw_curves): 1,2,4,6

# for old-style allcurves file
#awk '{print $1,$2,$3,$4;}' ${d}/allcurves/allcurves.${r} > temp.${r}
# for new-style allcurves file
awk '{print $1,$2,$4,$6;}' ${d}/allcurves/allcurves.${r} > temp.${r}

magma -b filename:=temp.${r} $HOME/ecdata/scripts/2adic.m
mv 2adic.${r} ${d}/2adic/
/bin/rm temp.${r}

