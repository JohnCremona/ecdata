#!/bin/bash

if [ $# != 2 ]; then
 echo "Usage: ./make_2adic.sh <range> <src_dir>"
 exit 1
else
    r=$1
    d=$2
fi

echo 'running 2-adic Magma script on '${r} in ${d}
awk '{print $10,$2,$3,$11;}' ${d}/curvedata/curvedata.${r} > temp.${r}
magma -b filename:=temp.${r} $HOME/ecdata/scripts/2adic.m
mv 2adic.${r} ${d}/2adic/
pushd ${HOME}/galrep
awk '{print $1":"$11;}' ${d}/curvedata/curvedata.${r} > temp.${r}
echo "ComputeQGaloisImages(\"temp.$r\", \"galrep.$r\"); quit;" > magma_temp.m
cat magma_temp.m
magma $HOME/galrep/nfgalrep.m magma_temp.m > /dev/null
mv galrep.${r} ${d}/galrep/
rm -f temp.${r} magma_temp.m
popd

