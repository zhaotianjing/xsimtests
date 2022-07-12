#! /bin/bash

snpALL=164

for snp in $( eval echo {1..$snpALL} )
do
    sbatch snp_level_pblup.sbatch $snp
done
