#!/bin/bash

if [ $# != 2 ]; then
 echo "Usage: ./make_2adic.sh <range> <src_dir>"
 exit 1
else
    r=$1
    d=$2
fi

echo 'running 2-adic Magma script on '${r} in ${d}

# Here we need the label in 3 parts then the a-invariant list, space-separated
# - in curvedata files this is fields 10,2,3,11
# - in allcurves files (as output by process_raw_curves): 1,2,4,6

awk '{print $1,$2,$4,$6;}' ${d}/allcurves/allcurves.${r} > temp.${r}
#awk '{print $10,$2,$3,$11;}' ${d}/curvedata/curvedata.${r} > temp.${r}
magma -b filename:=temp.${r} $HOME/ecdata/scripts/2adic.m
mv 2adic.${r} ${d}/2adic/

pushd ${HOME}/galrep

# Here we need the complete label then the a-invariant list, ":"-separated
# - in curvedata files this is fields 1,11
# - in allcurves files (as output by process_raw_curves): 1+2+4,6

awk '{print $1$2$4":"$6;}' ${d}/allcurves/allcurves.${r} > temp.${r}
#awk '{print $1":"$11;}' ${d}/curvedata/curvedata.${r} > temp.${r}
echo "ComputeQGaloisImages(\"temp.$r\", \"galrep.$r\"); quit;" > magma_temp.m
cat magma_temp.m
magma $HOME/galrep/nfgalrep.m magma_temp.m > /dev/null
mv galrep.${r} ${d}/galrep/
rm -f temp.${r} magma_temp.m
popd

