#!/bin/bash
npoints=$1

echo "Sourcing setup scripts"
source /root-6.12.06-build/bin/thisroot.sh

tar -zxvf pMSSM_McMC_v1.1.tar.gz
cd pMSSM_McMC

./install-hb.sh 
export HiggsBounds_DIR=${PWD}/packages/higgsbounds/build/

./install-hs.sh 

outdir=output_${npoints}
mkdir ${outdir}
python mcmc.py -m new -o ${outdir} -n ${npoints}
